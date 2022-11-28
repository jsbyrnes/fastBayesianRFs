function [ model ] = gibbs_iteration_prior(model, allWfs, TD_parameters, index)

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

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %change to A/dt
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind
        evtind      = 1:nevt;

        %field = 'A';
        if rand() > 0.5

            field = 'A';

        else

            field = 'dt';

        end
        type  = 'all';
        scale = TD_parameters.ABstd;%same search size for each

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to B/fastdir
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind
        evtind      = 1:nevt;
        
        %field = 'B';
        if rand() > 0.5

            field = 'B';
            scale = TD_parameters.ABstd;%same search size for each

        else

            field = 'fast_dir';
            scale = pi;

        end

        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to rotation
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind
        evtind      = 1:nevt;
        
        field = 'fast_dir_rotation';
        scale = TD_parameters.rotation_std;

        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to orientations
    if index <= nsta && isempty(field)

        staind = index;%applies directly
        evtind = 1:nevt;
        
        field = 'sta_or';
        scale = TD_parameters.sta_err(index);

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
    if isempty(field) %%%%%NEVER TRIGGERS

        evtind = 1:nevt;
        staind = 1:nsta;

        field = 'f';
        scale = TD_parameters.f_range(2);

        type  = 'error';%not really but this is where it is evaluated

    end

    %number of samples scale to however many times you update these
    starting_mpst = model.llh/model.T + model.lp;

    %a local or heavy tailed search
    
    if rand()>0.5

        pert = exprnd(1);

    else

        pert = exprnd(TD_parameters.step_std);

    end

    modeln                = model;
    modeln.(field)(index) = modeln.(field)(index) + scale*pert*(-1).^randi([0, 1]);
    modeln                = apply_update(modeln, field, TD_parameters.t);
    modeln                = evaluate(modeln, allWfs, TD_parameters, evtind, staind, type);
    lpst                  = modeln.llh/model.T + modeln.lp;

    if (lpst - starting_mpst) > log(rand())

        model = modeln;

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

