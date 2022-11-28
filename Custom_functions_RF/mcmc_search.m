function [mout, history] = mcmc_search(allWfs, TD_parameters, starting_model)

    alpha = 0.15;

    [nevt, nsta] = size(allWfs);

    if nargin == 2

         for k = 1:TD_parameters.nchains
    
            model(k) = make_starting_models_GA(TD_parameters, allWfs);
            model(k).vector_mean = model(k).vector;
            model(k).vector_std  = ones(size(model(k).vector));
    
        end

    else

         for k = 1:TD_parameters.nchains
    
            model(k)             = starting_model;            
            model(k).vector      = model(k).vector + 1e-1*randn(size(model(k).vector));
            model(k).vector_mean = model(k).vector;
            model(k).vector_std  = 1e-2*ones(size(model(k).vector));
    
            model(k) = devectorize_model(model(k), nevt, nsta, length(TD_parameters.t), TD_parameters);
            model(k) = evaluate(model(k), allWfs, TD_parameters, 1:nevt, 1:nsta, []);

         end

    end

    if TD_parameters.parallel
    
        if isempty(gcp('nocreate'))
    
            parpool('local');
    
        end
    
        disp('Distributing data to each worker before running the program')
        allWfs_const = parallel.pool.Constant(allWfs);
    
    end
    
    history(:, 1) = [ model(:).lpst ];
    mout          = [];
    
    for iter = 1:TD_parameters.iter
    
        if TD_parameters.parallel
    
            parfor j = 1:TD_parameters.nchains
            
                v0 = model(j).vector_mean;
                s0 = model(j).vector_std;
    
                for i = 1:TD_parameters.batch_size
    
                    model(j) = mcmc_iteration(model(j), allWfs, TD_parameters);
    
                end
    
                del                  = model(j).vector - v0;
                model(j).vector_mean = v0 + del*alpha;
                model(j).vector_std  = (1-alpha)*(s0 + alpha*del.^2);
    
            end
    
        else
    
            for j = 1:TD_parameters.nchains
            
                v0 = model(j).vector_mean;
                s0 = model(j).vector_std;
    
                for i = 1:TD_parameters.batch_size
    
                    model(j) = mcmc_iteration(model(j), allWfs, TD_parameters);
    
                end
    
                del                  = model(j).vector - v0;
                model(j).vector_mean = v0 + del*alpha;
                model(j).vector_std  = (1-alpha)*(s0 + alpha*del.^2);
        
            end
    
        end
    
        s0 = zeros(size(model(1).vector));

        for j = 1:TD_parameters.nchains
        
            s0 = s0 + log10(model(j).vector_std);
    
        end

        s0 = 10.^(s0/TD_parameters.nchains);

        for j = 1:TD_parameters.nchains
        
            model(j).vector_std = s0;
    
        end

        if TD_parameters.solver_printout
    
            disp([ 'At Iteration ' num2str(iter) ' with llh ' num2str(mean([ model.lpst ]))]);
            
            [~, n] = size(history);
            history(:, n+1) = [ model(:).lpst ];
    
        end
    
        if iter > TD_parameters.burnin && mod(iter, TD_parameters.batch_thin)
    
            msave = [];
    
            for k = 1:length(model)
                
                m = rmfield(model(k), 'Cinv');
                m = rmfield(m,        'irF');
                m = rmfield(m,        'irS');
                m = rmfield(m,        'N'  );
                m = rmfield(m,        'E'  );
    
                msave = [ msave; m ];
            
            end
    
            mout  = [ mout msave ];
            msave = [];
    
        end
    
    end

end

function model = mcmc_iteration(model, allWfs, TD_parameters)

    [nevt, nsta] = size(allWfs);
    model        = vectorize_model(model, TD_parameters);
    ind          = randperm(length(model.vector));
    modeln       = model;

    for k = 1:round(length(ind))

        modeln = coordinate_pert(modeln, allWfs, TD_parameters, ind(k), 1);

    end

    modeln = evaluate(modeln, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

    if (modeln.lpst - model.lpst) > log(rand())
       
        modeln            = vectorize_model(modeln, TD_parameters);
        %modeln.vector_del(ind) = modeln.vector_del(ind)*0.99 ...
        %    + 0.01*abs(modeln.vector(ind) - model.vector(ind));
        
        %del               = modeln.vector - model.vector_mean;
        model             = modeln;
        %model.vector_mean = model.vector_mean + del*alpha;
        %model.vector_std  = (1-alpha)*(model.vector_std + alpha*del.^2);

    else

        m1 = modeln;
        a1 = exp(modeln.lpst - model.lpst);
    
        m2 = model;

        sig2 = 0.2;

        q = 0;

        for k = 1:(length(ind))
    
            m2 = coordinate_pert(m2, allWfs, TD_parameters, ind(k), sig2);
            m2 = vectorize_model(m2, TD_parameters);

            q = q - ((m2.vector(ind(k)) - m1.vector(ind(k)))^2 - ...
            (m1.vector(ind(k)) - model.vector(ind(k)))^2)/(2*sig2^2);

        end
        
        m2 = evaluate(m2, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

        amid = min([1 exp(m1.lpst - m2.lpst) ]);

        q = exp(q);
    
        a = min([ 1 exp(m2.lpst - model.lpst)*q*(1 - amid)/(1 - a1)]);
    
        if rand() < a

            m2            = vectorize_model(m2, TD_parameters);

            %del               = m2.vector - model.vector_mean;
            model             = m2;
            %model.vector_mean = model.vector_mean + del*alpha;
            %model.vector_std  = (1-alpha)*(model.vector_std + alpha*del.^2);

        end
            
    end

end

function model = coordinate_pert(model, allWfs, TD_parameters, index, total_scale)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get the linear index of what will be changed
    [nevt, nsta] = size(allWfs);
    nwavelet     = numel(model.wavelet);%n in wavelets
    nsplit       = numel(model.A);
    scale  = 1;

    index0 = index;

    if nargin < 4
        index  = randi([1 model.nparam]);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %map the index to the structure fields
    evtind = [];
    staind = [];
    field  = [];

    %change to the wavelet
    if index <= nwavelet

        [evtind, ~] = ind2sub(size(model.wavelet), index);%check if this wraps right
        staind      = 1:nsta;
        
        field = 'wavelet';

    else

        index = index - nwavelet;

    end

    %change to polarization
    if index <= nevt && isempty(field)

        staind = 1:nsta; 
        evtind = index;%maps directly
        field  = 'polarization';
        scale  = TD_parameters.polarization_std;

    elseif isempty(field)

        index = index - nevt;

    end

    %%%%amp
    if index <= nsta*nevt && isempty(field)

        [evtind, staind] = ind2sub(size(allWfs), index);
        
        field = 'amp';
        scale = 1;

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %%%%delay
    if index <= nsta*nevt && isempty(field)

        [evtind, staind] = ind2sub(size(allWfs), index);
        
        field = 'shift';
        scale = TD_parameters.Delaystd;

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %%%%tS
    if index <= nsta*nevt && isempty(field)

        [evtind, staind] = ind2sub(size(allWfs), index);
        
        field = 'dtS';
        scale = TD_parameters.dtS(2);

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %change to A
    if index <= nsplit && isempty(field)

        if TD_parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end

        evtind      = 1:nevt;

        field = 'A';
        scale = TD_parameters.ABstd;%same search size for each

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to B
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        if TD_parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end

        evtind      = 1:nevt;
            
        field = 'B';
        scale = TD_parameters.ABstd;%same search size for each

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to rotation
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        [~, lind] = ind2sub(size(model.A), index);
        if TD_parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end
        evtind = 1:nevt;
        field  = 'fast_dir_rotation';
        scale  = TD_parameters.rotation_std;

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

    elseif isempty(field)

        index = index - nsta*nevt;

    end
        
    %change to r
    if index == 1 && isempty(field)

        index = 1;%scalar
        
        evtind = 1:nevt;
        staind = 1:nsta;
        scale = TD_parameters.r_range(2);
        field = 'r';

    elseif isempty(field)

        index = index - 1;

    end

    %change to f
    if index == 1 && isempty(field)

        evtind = 1:nevt;
        staind = 1:nsta;
        field  = 'f';
        scale  = TD_parameters.f_range(2);

    elseif isempty(field)

        index = index - 1;

    end

    scale = total_scale*scale;%used for delayed rejection

    %model.(field)(index) = model.(field)(index) + TD_parameters.pert*scale*randn();
    model.(field)(index) = model.(field)(index) + min([ sqrt(model.vector_std(index0)) TD_parameters.pert ])*scale*randn();
    %model.(field)(index) = model.(field)(index) + sqrt(model.vector_std(index0))*scale*randn();
    model = apply_update(model, field, TD_parameters.t);

end