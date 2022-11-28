function [ mout, history ] = fit_waveforms_old(allWfs, TD_parameters, synmodel)
%synmodel only used in debugging

%iwth a bunch of extra shit that I tried during development

    [nevt, nsta] = size(allWfs);

    if strcmp(TD_parameters.solver, 'bfgs')

        history = 0;

        disp('-> Building a starting model')
        model = make_starting_models_GA(TD_parameters, allWfs);

        model = vectorize_model(model, TD_parameters);

        if TD_parameters.solver_printout

            options                  = optimoptions('fminunc','Display','iter','Algorithm','quasi-newton',...
                'MaxIterations', 1e9, 'OptimalityTolerance', 1e-4, ...
                'MaxFunctionEvaluations', 1e9, 'SpecifyObjectiveGradient', true);
         
        else

            options                  = optimoptions('fminunc','Display','none','Algorithm','quasi-newton',...
                'MaxIterations', 1e9, 'OptimalityTolerance', 1e-4, ...
                'MaxFunctionEvaluations', 1e9, 'SpecifyObjectiveGradient', true);
            
        end

        %first attempt at fitting
        disp('-> Fitting a first model')

        f                        = @(x) fitwaveform_wrapper(x, model, allWfs, TD_parameters, -1);
        [x, ~, ~, fmin_out, del] = fminunc(f, model.vector, options);

        model.vector   = x;
        model          = devectorize_model(model, nevt, nsta, length(TD_parameters.t), TD_parameters);
        model          = evaluate(model, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
        model.del      = del;
        model.fmin_out = fmin_out;

        disp([ '--> First solution at posterior density of ' num2str(model.lpst) ]);

        iter         = 1;
        rounds       = 0;
        reset_thresh = 0.1;%not important to fine tune, threshold for when you don't count as an improvement

        modelr = model;%for the model in this round

        while rounds < TD_parameters.reset_rounds
            
            if iter == 0

                modelr = make_starting_models_GA(TD_parameters, allWfs);
                modelr = vectorize_model(modelr, TD_parameters);
            
                iter = 1;

            end
        
            sig           = 10^(diff(TD_parameters.reset_size)*(iter/TD_parameters.reset_iter)...
                + TD_parameters.reset_size(1));
            modeln        = modelr;
            modeln.vector = modelr.vector + sig*randn(size(modelr.vector));
            modeln        = devectorize_model(modeln, nevt, nsta, length(TD_parameters.t), TD_parameters);

            modeln        = evaluate(modeln, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
            disp([' -> Attempting a trial starting at posterior density of ' num2str(modeln.lpst) ]);
            
            f                        = @(x) fitwaveform_wrapper(x, modeln, allWfs, TD_parameters, -1);
            [x, ~, ~, fmin_out, del] = fminunc(f, modeln.vector, options);

            modeln.vector   = x;
            modeln          = devectorize_model(modeln, nevt, nsta, length(TD_parameters.t), TD_parameters);
            modeln          = evaluate(modeln, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

            if (modeln.lpst - modelr.lpst) > reset_thresh

                disp([' -> Better trial solution found with posterior density of ' num2str(modeln.lpst) ]);
                modelr          = modeln;
                modelr.del      = del;
                modelr.fmin_out = fmin_out;

                iter = max([1 (iter - 1)]);

            elseif reset_thresh > (modeln.lpst - modelr.lpst) && (modeln.lpst - modelr.lpst) >= 0

                disp( '    Solution found within search tolerence');
                iter = iter + 1;
                disp(['    Counter raised to ' num2str(iter) ' of ' num2str(TD_parameters.reset_iter) ])

            else

                disp(['    Test solution rejected at posterior density of ' num2str(modeln.lpst) '']);                
                iter = iter + 1;
                disp(['    Counter raised to ' num2str(iter) ' of ' num2str(TD_parameters.reset_iter) ])

            end

            if iter == TD_parameters.reset_iter

                if (modelr.lpst - model.lpst) > reset_thresh

                    disp(['--> Better solution over a round found with posterior density of ' num2str(modelr.lpst) ]);
                    disp( '--> **Round counter set to zero**');
                    
                    model = modelr;
    
                    rounds = 0;

                else

                    rounds = rounds + 1;

                    disp('    ')
                    disp([ '--> Round counter raised to ' num2str(rounds)])

                end

                modelr = model;%and then reset it
                iter   = 0;%start over

            end

        end

        %uncertainities!
        if strcmp(TD_parameters.error_type, 'hessian')

            %total bullshit!
            disp('-> Getting the Hessian for uncertainities')
            model   = fill_hessian(model, allWfs, TD_parameters);
            model.C = inv(-model.H);%covarience matrix
            mout    = model;

        elseif strcmp(TD_parameters.error_type, 'hmc')

            %seems to work the best
            f  = @(x) fitwaveform_wrapper(x, model, allWfs, TD_parameters, 1);

            smp = hmcSampler(f, model.vector, 'StepSize', TD_parameters.pert, ...
                'NumSteps', 5, 'MassVectorTuningMethod', 'hessian', 'StepSize', 1e-1);

            c = [];

            disp('-> Starting HMC search for errors');
            if TD_parameters.solver_printout

                smp          = tuneSampler(smp, 'VerbosityLevel', 2, 'NumPrint',1, 'NumStepSizeTuningIterations', 25, 'NumStepsLimit', 5);
                smp.NumSteps = 5;%often diverges during tuning despite numstepslimits flag

                for k = 1:TD_parameters.nchains

                    tmp = drawSamples(smp, 'Burnin', TD_parameters.burnin, 'NumSamples', TD_parameters.batch_size, ...
                        'ThinSize', TD_parameters.batch_thin ...
                    , 'Verbosity', 1, 'NumPrint', 1, 'StartPoint', model.vector);

                    c = [ c; tmp ];

                end

            else

                smp          = tuneSampler(smp, 'NumStepSizeTuningIterations', 25, 'NumStepsLimit', 5);
                smp.NumSteps = 5;%often diverges during tuning despite numstepslimits flag

                for k = 1:TD_parameters.nchains

                    disp([ '     Chain #' num2str(k) ' of ' num2str(TD_parameters.nchains) ]);

                    tmp = drawSamples(smp, 'Burnin', TD_parameters.burnin, 'NumSamples', TD_parameters.batch_size, ...
                        'ThinSize', TD_parameters.batch_thin, 'StartPoint', model.vector);

                    c = [ c; tmp ];

                end

            end

            disp('   Building models from HMC search');
            for k = 1:length(c)

                m        = model;
                m.vector = c(k, :);
                m        = devectorize_model(m, nevt, nsta, length(TD_parameters.t), TD_parameters);
                mout(k)  = evaluate(m, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

            end
                    
        end

    elseif strcmp(TD_parameters.solver, 'acd') %often converges prematurally

        history = 0;
        
        disp('-> Building a starting model')
        model = make_starting_models_GA(TD_parameters, allWfs);
        model = vectorize_model(model, TD_parameters);

        if TD_parameters.parallel
    
            xmean = ACD( @(XV,mu) ParForParallelWrapper(f, XV, mu, 0 ), model.vector,5,1e-6,[],[],[],[],[],[],1,...
                2,2,2,TD_parameters.printout);
        
        else
    
            xmean = ACD( @(XV,mu) SerialWrapper( f, XV, mu, 0 ), model.vector,5,1e-6,[],[],[],[],[],[],2,...
                2,2,2,TD_parameters.printout);
    
        end
    
%         [xmean, ~, ~] = ACD(f,length(model.vector),-25,25,1e12,-2e4,1,model.vector);
% 
        model.vector = xmean;
        model        = devectorize_model(model, nevt, nsta, length(TD_parameters.t), TD_parameters);
        model        = evaluate(model, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

        mout = model;

    elseif strcmp(TD_parameters.solver, 'mcmc')

        [ mout, history ] = mcmc_search(allWfs, TD_parameters);

    elseif strcmp(TD_parameters.solver, 'hmc')

        %may not converge to a reasonable
        %optimum!

        c = [];

        mout = [];

        disp('-> Starting HMC search from random models');

        if TD_parameters.solver_printout

            for k = 1:TD_parameters.nchains

                model = make_starting_models_GA(TD_parameters, allWfs);
                f     = @(x) fitwaveform_wrapper(x, model, allWfs, TD_parameters, 1);
                smp   = hmcSampler(f, model.vector, 'StepSize', 0.005, 'NumSteps', 5);
                smp   = tuneSampler(smp, 'VerbosityLevel', 2, 'NumPrint',1, 'NumStepSizeTuningIterations', 25);
                tmp   = drawSamples(smp, 'Burnin', TD_parameters.burnin, 'NumSamples', TD_parameters.batch_size, ...
                    'ThinSize', TD_parameters.batch_thin ...
                , 'Verbosity', 1, 'NumPrint',1);

                for j = 1:length(tmp)

                    m        = model;
                    m.vector = tmp(j, :);

                    mout = [ mout evaluate(devectorize_model(m, nevt, nsta, length(TD_parameters.t), TD_parameters), allWfs, TD_parameters, 1:nevt, 1:nsta, []) ];

                end

            end

        else

            for k = 1:TD_parameters.nchains

                disp([ '     Chain #' num2str(k) ' of ' num2str(TD_parameters.nchains) ]);

                model = make_starting_models_GA(TD_parameters, allWfs);
                f     = @(x) fitwaveform_wrapper(x, model, allWfs, TD_parameters, 1);
                smp   = hmcSampler(f, model.vector, 'StepSize', 0.005, 'NumSteps', 5);
                smp   = tuneSampler(smp, 'NumStepSizeTuningIterations', 25);
                tmp   = drawSamples(smp, 'Burnin', TD_parameters.burnin, 'NumSamples', TD_parameters.batch_size, ...
                    'ThinSize', TD_parameters.batch_thin);

                for j = 1:length(tmp)

                    m        = model;
                    m.vector = tmp(j, :);
                    
                    mout = [ mout evaluate(devectorize_model(m, nevt, nsta, length(TD_parameters.t), TD_parameters), allWfs, TD_parameters, 1:nevt, 1:nsta, []) ];

                end

            end

        end

        model.vector_error = std(c);
        mout               = model;

%     elseif strcmp(TD_parameters.solver, 'gibbs')
% 
%         disp('-> Building a starting model')
%         model = make_starting_models_GA(TD_parameters, allWfs);
% 
%         model = vectorize_model(model, TD_parameters);
% 
%         pert  = 1e9*ones(model.nparam, 1);
%         niter = 0;
% 
%         while niter < 1000 && max(abs(pert)) > 0.001
%         
%             inds  = randperm(model.nparam);
%             for k = 1:length(inds)
%                [model, pert(k)] = coordinate_descent(model, allWfs, TD_parameters, inds(k));
%             end        
%         
%             niter = niter + 1;
%     
%             disp([ 'Iteration iteration ' num2str(niter) ' with max ' num2str(max(abs(pert))) ...
%                 ' and llh ' num2str(model.llh)]);
%     
%         end

    end

end

% function [model, pert] = coordinate_descent(model, allWfs, TD_parameters, index)
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %get the linear index of what will be changed
%     [nevt, nsta] = size(allWfs);
%     nwavelet     = numel(model.wavelet);%n in wavelets
%     nsplit       = numel(model.A);
%     scale        = 1;
% 
%     index0 = index;%for debugging
% 
%     if nargin < 4
%         index  = randi([1 model.nparam]);
%     end
% 
%     [index, field, scale, evtind, staind] = map_index(index0, model, allWfs, TD_parameters);
% 
%     %%%%%%
%     %based on the inverse problem for alignment
%     pert  = 0;
%     diff  = Inf;
%     niter = 0;
% 
%     lpst = model.llh + model.lp;
% 
%     while diff > h
% 
%         lpst0 = make_new(model, field, index, pert,  allWfs, TD_parameters, evtind, staind);
%         lpstp = make_new(model, field, index, pert + h,  allWfs, TD_parameters, evtind, staind);
%         lpstn = make_new(model, field, index, pert - h, allWfs, TD_parameters, evtind, staind);
% 
%         d1 = (lpstp - lpstn)/(2*h);
%         d2 = (lpstp + lpstn - 2*lpst0)/(h^2);
% 
%         step  = d1/d2;
% 
%         if isnan(step) && d2 == 0
% 
%             break
% 
%         end
% 
%         if d2 < 0 && abs(step) <= 1 %blows up when second derivative is small
% 
%             pertn = pert - step;
%             diff  = abs(pertn - pert);
%             pert  = pertn;
% 
%         elseif (niter ~= 10) && (d2 > 0 || abs(step) >= 1) %newton's method won't work! psuedo grid search to reset the parameter
% 
%             lpst0   = make_new(model, field, index, 0,  allWfs, TD_parameters, evtind, staind);
%             alpha  = logspace(-2, 0, 10);
%             maxlpst = -Inf;
% 
%             while maxlpst < lpst0
% 
%                 %line search gradient descent
%                 step_v     = sign(d1)*alpha;
%                 pert_v     = pert + step_v;
%                 lpst_search = [];
%                 for q = 1:length(pert_v)
%                     lpst_search(q) = make_new(model, field, index, pert_v(q),  allWfs, TD_parameters, evtind, staind);
%                 end
% 
%                 [maxlpst, ii] = max(lpst_search);
% 
%                 if maxlpst > lpst0
% 
%                     diff  = abs(pert - pert_v(ii));
%                     pert = pert_v(ii);%otherwise, marches into shit territory
% 
%                 end
% 
%                 alpha = alpha/5;
% 
%                 if any(alpha < 1e-12)
% 
%                     %pathological, so give up and grid
%                     %search
%                     %keyboard
%                     x = (-(3*scale):h:(3*scale));
%                     lpst_search = [];
%                     for q = 1:length(x)
%                         lpst_search(q) = make_new(model, field, index, x(q),  allWfs, TD_parameters, evtind, staind);
%                     end
% 
%                     [~, id] = max(lpst_search);
%                     pert    = x(id);
%                     break
% 
%                 end
% 
%             end
% 
%         end
% 
%         if niter == 10
% 
%             x = -(2*scale):(scale/10):(2*scale);
%             lpst_search = [];
%             for q = 1:length(x)
%                 lpst_search(q) = make_new(model, field, index, x(q),  allWfs, TD_parameters, evtind, staind);
%             end
% 
%             [~, id] = max(lpst_search);
%             pert = x(id);
% 
%         end
% 
%         niter  = niter + 1;
% 
%         if niter > 50
% 
%             %pathological, so give up and grid
%             %search
% 
%             disp(['Grid search required for ' field])
% 
%             x = (-(3*scale):(h*10):(3*scale));
%             lpst_search = [];
%             for q = 1:length(x)
%                 lpst_search(q) = make_new(model, field, index, x(q),  allWfs, TD_parameters, evtind, staind);
%             end
% 
%             [~, id] = max(lpst_search);
%             pert    = x(id);
%             break
% 
%         end
% 
%     end
% 
%     [~, model] = make_new(model, field, index, pert,  allWfs, TD_parameters, evtind, staind);
% 
%     %normalize pert to scale
%     pert = pert/scale;
%         
% end
% 
% function [lpt, model] = make_new(model, field, index, pert, allWfs, TD_parameters, evtind, staind)
% 
%     model.(field)(index) = model.(field)(index) + pert;
% 
%     if field == 'f'
% 
%         model.(field)(index) = min([model.(field)(index) TD_parameters.sample_rate]);
% 
%     end
% 
%     if strcmp(field, 'sig')
% 
%         model.(field)(index) = min([model.(field)(index) 1]);
% 
%     end
% 
%     model = apply_update(model, field, TD_parameters.t);
%     model = evaluate(model, allWfs, TD_parameters, evtind, staind, []);
% 
%     lpt = model.lpst;
% 
% end
% 
%     exitflag = -1e9;
% 
%     while (exitflag ~= 1)
% 
%         %premature exit often caused by huge gradients in f/r, so optimize
%         %when this happens and try again. Ok to do several iterations
%         model = center_wavelets(model, allWfs, TD_parameters);
% 
%         inds  = randperm(model.nparam);
%         for k = 1:length(inds)
%            [model, ~] = coordinate_descent(model, allWfs, TD_parameters, inds(k));
%         end        
%         model = newton_C(model, allWfs, TD_parameters);
%         model = vectorize_model(model, TD_parameters);
%         try
% 
%             [x,~,exitflag,~] = fminunc(f, model.vector, options);
%             model.vector = x;
%             model        = devectorize_model(model, nevt, nsta, length(TD_parameters.t), TD_parameters);
%             model        = evaluate(model, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
% 
%         catch
% 
%             exitflag = -1e9;%errors when rounnding means its not a descent direction
% 
%         end
% 
%         disp([ num2str(model.fast_dir*180/pi), ';', num2str(model.fast_dir_rotation*180/pi) ';' num2str(model.dt) ';' num2str(model.llh)])
% 
%     end
% 
% %     options                  = optimset('Display', 'iter', 'MaxIter', 10000);
% % 
% %     [x,fval,exitflag,output] = fminsearch(f, model.vector, options);
% 

%     [xmean, fcurrent, neval] = ACD(f,length(model.vector),-5,5,1e12,-2e4,length(model.vector)/10,model.vector);
% 
%     model.vector = xmean;
%     model        = devectorize_model(model, nevt, nsta, length(TD_parameters.t), TD_parameters);
%     model        = evaluate(model, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

%         elseif strcmp(TD_parameters.error_type, 'mcmc')
% 
%             %works, but super inefficient!
%             mout = mcmc_search(allWfs, TD_parameters, model);


