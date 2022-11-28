function [ model ] = gibbs_iteration(model, allWfs, TD_parameters, index)

    mtmp = model;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get the linear index of what will be changed
    [nevt, nsta] = size(allWfs);
    nwavelet     = numel(model.wavelet);%n in wavelets
    nsplit       = numel(model.A);%number of splitting parameters. Could be clustered
    
    if nargin < 4
        index  = randi([1 model.nparam]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %map the index to the structure fields
    evtind = [];
    staind = [];
    type   = [];
    field  = [];
    scale  = 1;
    joint  = false;

    upper_bound = Inf;%some parameters need hard bounds

    %change to the wavelet
    if index <= nwavelet

        [evtind, ~] = ind2sub(size(model.wavelet), index);
        staind      = 1:nsta;
        
        field = 'wavelet';
        type  = 'all';

    else

        index = index - nwavelet;

    end

    %change to polarization
    if index <= nevt && isempty(field)

        staind = 1:nsta; 
        evtind = index;%maps directly
        scale  = TD_parameters.polarization_std;
        field  = 'polarization';
        type   = 'all';

    elseif isempty(field)

        index = index - nevt;

    end

    %change to dtS
    if index <= nsta*nevt && isempty(field)

        [evtind, staind] = ind2sub(size(allWfs), index);
        
        field = 'dtS';
        type  = 'all';
        scale = TD_parameters.dtS(2);
        upper_bound = log(5);

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %change to A/dt
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        [~, lind] = ind2sub(size(model.A), index);%dont need layer ind

        if TD_parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end

        evtind      = 1:nevt;

        %field = 'A';

        if TD_parameters.prior_information.use

            field       = 'dt';
            upper_bound = 5;
            scale       = TD_parameters.prior_information.dt_std(lind);

        else

            if rand() > 0.5
    
                field = 'A';
    
            else
    
                field = 'dt';
                upper_bound = 5;
    
            end
            scale = TD_parameters.ABstd;%same search size for each

        end

        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to B/fastdir
    if index <= nsplit && isempty(field)

        [~, lind] = ind2sub(size(model.A), index);%dont need layer ind

        %A works for size on all splitting parameters
        if TD_parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end
        evtind      = 1:nevt;
        
        if TD_parameters.prior_information.use

            field       = 'fast_dir';
            scale       = TD_parameters.prior_information.phi_std(lind);

        else

            if rand() > 0.5
    
                field = 'B';
                scale = TD_parameters.ABstd;%same search size for each
    
            else
    
                field = 'fast_dir';
                scale = pi;
    
            end

        end

        %field = 'B';
        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to rotation
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        [~, lind] = ind2sub(size(model.A), index);%dont need layer ind
        if TD_parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end
        evtind = 1:nevt;
        field  = 'fast_dir_rotation';

        if TD_parameters.prior_information.use

            scale       = TD_parameters.prior_information.rot_std(lind);

        else

            scale = TD_parameters.rotation_std;

        end

        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to orientations
    if index <= nsta && isempty(field)

        staind = index;%applies directly
        evtind = 1:nevt;
        
        field = 'sta_or';
        scale = min([ TD_parameters.sta_err(index), pi/2]);

        type  = 'all'; 

    elseif isempty(field)

        index = index - nsta;

    end
    
    %change to sigma
    if index <= nsta*nevt && isempty(field)

        [evtind, staind] = ind2sub(size(allWfs), index);
        
        field = 'sig';
        scale = TD_parameters.sig_range(2);

        type  = 'error';%not really but this is where it is evaluated
        upper_bound = 1;

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %change to r
    if index == 1 && isempty(field)

        index = 1;%scalar
        
        evtind = 1:nevt;
        staind = 1:nsta;

        field = 'r';
        scale = TD_parameters.r_range(2);

        type  = 'error';%not really but this is where it is evaluated

    elseif isempty(field)

        index = index - 1;

    end

    %change to f
    if index == 1

        evtind = 1:nevt;
        staind = 1:nsta;

        field       = 'f';
        scale       = TD_parameters.f_range(2);
        upper_bound = TD_parameters.low_pass;

        type  = 'error';%not really but this is where it is evaluated

    elseif isempty(field)

        index = index - 1;

    end

    %number of samples scale to however many times you update these
    starting_mpst = model.llh/model.T + model.lp;

    %mix of a local and heavy tailed search
    
    pert = scale*[ exprnd(1, round(TD_parameters.line_iterations), 1)...
        .*(-1).^randi([0, 1], round(TD_parameters.line_iterations), 1);
        exprnd(TD_parameters.step_std, round(TD_parameters.line_iterations), 1)...
        .*(-1).^randi([0, 1], round(TD_parameters.line_iterations), 1) ];

    for kk = 1:length(pert)

        modeln                = model;
        modeln.(field)(index) = modeln.(field)(index) + pert(kk);

        if modeln.(field)(index) >= upper_bound

            lpst(kk) = -Inf;

        else

            modeln                = apply_update(modeln, field, TD_parameters.t);
            modeln                = evaluate(modeln, allWfs, TD_parameters, evtind, staind, type);
            lpst(kk)              = modeln.llh/model.T + modeln.lp;

        end

    end

    %this is useful when debugging
    %plot(pert + model.(field)(index), lpst - starting_mpst, 'k.'), hold on, ...
    %plot(model.(field)(index), model.llh + model.lp - starting_mpst, 'k*')

%     x = -10:0.1:10;
%     y = a*(x - pert(i) - model.(field)(index)).^2 + b*(x - pert(i) - model.(field)(index)) + (mlpst - starting_mpst);
%     hold on
%     plot(x, y)

    [mlpst, i] = max(lpst);

%         %calculate relative probablities and make a random selection
     pert  = [ pert; 0 ];%might stay where you are
     lpst  = [ lpst starting_mpst];

    %three point parabola

    h = scale*TD_parameters.step_std;

    modelp                = model;
    modelp.(field)(index) = model.(field)(index) + pert(i) + h;
    modelp                = apply_update(modelp, field, TD_parameters.t);
    modelp                = evaluate(modelp, allWfs, TD_parameters, evtind, staind, type);
    
    modeln                = model;
    modeln.(field)(index) = model.(field)(index) + pert(i) - h;
    modeln                = apply_update(modeln, field, TD_parameters.t);
    modeln                = evaluate(modeln, allWfs, TD_parameters, evtind, staind, type);

    if ~((mlpst > ((modelp.llh + modelp.lp))) && (mlpst > ((modeln.llh + modeln.lp))))

        %increase the step size to sample the parabola
        h = scale*TD_parameters.step_std*3;
    
        modelp                = model;
        modelp.(field)(index) = model.(field)(index) + pert(i) + h;
        modelp                = apply_update(modelp, field, TD_parameters.t);
        modelp                = evaluate(modelp, allWfs, TD_parameters, evtind, staind, type);
        
        modeln                = model;
        modeln.(field)(index) = model.(field)(index) + pert(i) - h;
        modeln                = apply_update(modeln, field, TD_parameters.t);
        modeln                = evaluate(modeln, allWfs, TD_parameters, evtind, staind, type);

    end

    if ~((mlpst > ((modelp.llh + modelp.lp))) && (mlpst > ((modeln.llh + modeln.lp)))) %still a no?

        if mlpst > starting_mpst %just keep best. Must be shit location. 
    
            mtmp.(field)(index) = pert(i) + model.(field)(index);
            mtmp                = apply_update(mtmp, field, TD_parameters.t);

        end

    else

        yp = modelp.llh + modelp.lp - mlpst;
        yn = modeln.llh + modeln.lp - mlpst;

        b = (yp - yn)/(2*h);
        a = (yp + yn)/(2*h^2);

        optimum = -b/(2*a) + pert(i) + model.(field)(index);

        if TD_parameters.optimize

            mtmp.(field)(index) = optimum;

        else

            %quadratic formula at a random change to in log posterior
            %newpert = (-b + sqrt(b^2 + 4*a*log(rand()))*(-1)^randi([0 1]))/(2*a);%this is wrong

            %%%approximation alert - don't extrapolate with this
            %%%polynomial.             
            %newpert = min([ max([ min(pert) newpert ]) max(pert) ]);

            %convert the parabola to a distribution draw a random number
            %from it by taking the cdf from 0 to 1. Higher order polys
            %could be asymmetric.
            x    = -10:0.01:10;%-10 log units is ~1e-5 probablity. reasonable clip
            dist = exp(a*x.^2 + b*x);
            dist(isinf(dist)) = max(dist(~isinf(dist)));
            cdf  = cumsum(dist)/sum(dist) + cumsum(1e-10*ones(size(dist)));

            newpert = interp1(cdf, x, rand());
            mtmp.(field)(index) = pert(i) + model.(field)(index) + newpert;            

        end

        mtmp                = apply_update(mtmp, field, TD_parameters.t);

    end

    mtmp = evaluate(mtmp, allWfs, TD_parameters, evtind, staind, type);%burn a run
    
    if (mtmp.llh + mtmp.lp - starting_mpst) > 0 && TD_parameters.optimize

        model = mtmp;

    elseif (mtmp.llh + mtmp.lp - starting_mpst) > -10 && ~TD_parameters.optimize

        model = mtmp;

    else

        if mtmp.llh > 1000

            ffff=1;%so I can trigger it if I want
    
        end

        %failed parabola - optimize
        if mlpst > starting_mpst %just keep best

            model.(field)(index) = pert(i) + model.(field)(index);
            model                = apply_update(model, field, TD_parameters.t);
            model                = evaluate(model, allWfs, TD_parameters, evtind, staind, type);%burn a run

        end

    end

end

function model = apply_update(model, field, t)

    if strcmp(field, 'A') || strcmp(field, 'B')

        model = make_dtphi(model);

    elseif strcmp(field, 'dt') || strcmp(field, 'fast_dir')
    
        model = make_AB(model);
        model = make_dtphi(model);%does the wrapping correctly for fast_dir

    elseif strcmp(field, 'r') || strcmp(field, 'f')

        model = build_C(model, t);

    end

end


%%%%%%
%Old selection proceedure. Generated numerical errors of unknown cause.
%Only one parameters gets perturbed, but many parameters gets changes by a
%small amount. 
%     modeln               = model;
%     modeln.vector(index) = modeln.vector(index) + 0.5;
%     modeln               = devectorize_model(modeln, nevt, nsta, n, nl, TD_parameters);
% 
%     evtind = [];
%     staind = [];
%     type   = [];
% 
%     if index > (length(model.vector) - 2)
% 
%         type = 'error';
%     
%         evtind = 1:nevt;
%         staind = 1:nsta;
% 
%         fordebugging = 'covarience';
% 
%     elseif any(modeln.sig(:) - model.sig(:))
% 
%         type = 'error';
% 
%         [~, k]           = max(abs(modeln.sig(:) - model.sig(:)));%now idea why but two variables can change
%         [evtind, staind] = ind2sub([nevt, nsta], k);
% 
%         fordebugging = 'sig';
% 
%     elseif any(abs(circ_dist(modeln.sta_or,model.sta_or))>1e-10)
% 
%         type = 'error';
% 
%         [~, staind] = max(abs(modeln.sta_or - model.sta_or));
%         evtind = 1:nevt;
% 
%         fordebugging = 'or';
% 
%     elseif any(abs(circ_dist(modeln.polarization,model.polarization))>1e-10)
% 
%         type = 'all';
% 
%         [~, evtind] = max(abs(modeln.polarization - model.polarization));%rounding errors occur
%         staind = 1:nsta;
% 
%         fordebugging = 'polarization';
% 
%     elseif (any(modeln.A(:) - model.A(:))...
%             || any(modeln.B(:) - model.B(:))...
%             || any(modeln.fast_dir_rotation(:) - model.fast_dir_rotation(:))) ...
%             && model.vector_style == 0 %AB style
%     
%         type = 'all';
% 
%         if TD_parameters.cluster
% 
%             staind = 1:nsta;
% 
%         else
% 
%             staind = [ find((sum(modeln.A,2) - sum(model.A,2) )), ...
%                 find((sum(modeln.B,2) - sum(model.B,2) )), ...
%                 find((sum(modeln.fast_dir_rotation,2) - sum(model.fast_dir_rotation,2) ))];
% 
%         end
% 
%         evtind = 1:nevt;
% 
%         fordebugging = 'AB';
% 
%     elseif (any( abs(modeln.dt(:) - model.dt(:))>1e-10)...
%             || any( abs(modeln.fast_dir(:) - model.fast_dir(:))> 1e-10)...
%             || any(modeln.fast_dir_rotation(:) - model.fast_dir_rotation(:))) ...
%             && model.vector_style == 1 %spherical style. Note rounding error correction. For both log and atan
%     
%         type = 'all';
% 
%         if TD_parameters.cluster
% 
%             staind = 1:nsta;
% 
%         else
% 
%             staind = [ find(abs(sum(modeln.dt,2) - sum(model.dt,2) )>1e-10), ...
%                 find(abs(sum(modeln.fast_dir,2) - sum(model.fast_dir,2) ) > 1e-10), ...
%                 find((sum(modeln.fast_dir_rotation,2) - sum(model.fast_dir_rotation,2) ))];
% 
%         end
% 
%         evtind = 1:nevt;
% 
%         fordebugging = 'spherical';
% 
%     elseif any(modeln.wavelet(:) - model.wavelet(:))
%         
%         type = 'source_fit';
% 
%         [~, evtind] = max(sum(modeln.wavelet - model.wavelet, 2));
%         staind = 1:nsta;
%     
%         fordebugging = 'wavelet';
% 
%     elseif any(modeln.dtS(:) - model.dtS(:))
% 
%         type = 'source_fit';
% 
%         [~, k] = max(modeln.dtS(:) - model.dtS(:));
%         [evtind, staind] = ind2sub([nevt, nsta], k);
% 
%         fordebugging = 'sourcefitting';
% 
%     end

