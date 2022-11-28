function [ model, step_size, TD_parameters ] = gradient_ascent(model, allWfs, TD_parameters, step_size)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get the linear index of what will be changed
    [nevt, nsta] = size(allWfs);
    nwavelet     = numel(model.wavelet);%n in wavelets
    nsplit       = numel(model.A);%number of splitting parameters. Could be clustered

    starting_logp = model.llh + model.lp;
    step_v        = exp( log(step_size) + TD_parameters.step_rescaling*randn(TD_parameters.line_iterations,1) );

    for k = 1:TD_parameters.line_iterations

        modeln = model;

        step     = step_v(k)*model.del;

        %change to the wavelet
        modeln.wavelet   = model.wavelet + reshape(step(1:nwavelet), size(model.wavelet));
        step(1:nwavelet) = [];
        
        %change to polarization
        modeln.polarization = model.polarization + step(1:nevt);
        step(1:nevt)        = [];

        %change to dtS    
        modeln.dtS        = model.dtS + reshape(step(1:nsta*nevt), size(model.dtS));
        step(1:nsta*nevt) = [];

%         %change to A/dt
%         modeln.A       = model.A + reshape(step(1:nsplit), size(model.A));
%         step(1:nsplit) = [];
% 
%         %change to B/fastdir
%         modeln.B       = model.B + reshape(step(1:nsplit), size(model.A));
%         step(1:nsplit) = [];
%         modeln         = apply_update(modeln, 'A');
% 
%         %change to rotation
%         modeln.fast_dir_rotation = model.fast_dir_rotation + reshape(step(1:nsplit), size(model.A));
%         step(1:nsplit)           = [];
    
        %change to orientations
        modeln.sta_or = model.sta_or + reshape(step(1:nsta), size(model.sta_or));
        step(1:nsta)  = [];

        %change to sigma
        modeln.sig        = model.sig + reshape(step(1:nsta*nevt), size(model.sig));
        step(1:nsta*nevt) = [];
    
        %change to r
        modeln.r = model.r + step(1);
        step(1)  = [];
    
        %change to f
        modeln.f = log(min( [ exp(model.f + step(1)) TD_parameters.sample_rate]));
        modeln   = apply_update(modeln, 'r', TD_parameters.t);

        mset(k) = evaluate(modeln, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

    end

    llh = [mset(:).llh];
    lp  = [mset(:).lp];

    [mp, i] = max(llh + lp);

    if step_v(i) == min(step_v) 

        step_size = step_v(i)*(1 - 0.2*rand());
        TD_parameters.step_rescaling = TD_parameters.step_rescaling*(1 + 0.2*rand());

    elseif step_v(i) == max(step_v)

        step_size = step_v(i)*(1 + 0.2*rand());
        TD_parameters.step_rescaling = TD_parameters.step_rescaling*(1 + 0.2*rand());

    else

        TD_parameters.step_rescaling = TD_parameters.step_rescaling/(1 + 0.2*rand());
        
    end

    if (mp - starting_logp) > 0

        model     = mset(i);

    end    

end
