function mout = fit_waveforms_IDRF(allWfs, Parameters)
%synmodel only used in debugging

    tol = 0.1;%doesn't need to fine tuned

    [nevt, nsta] = size(allWfs);

    disp('-> Building a starting model')
    model = make_starting_models(Parameters, allWfs);

    if Parameters.solver_printout

%         options                  = optimoptions('fminunc','Display','iter','Algorithm','quasi-newton',...
%             'MaxIterations', 1e9, 'OptimalityTolerance', 1e-2, ...
%             'MaxFunctionEvaluations', 1e9, 'UseParallel', Parameters.parallel);
        options                  = optimoptions('fminunc','Display','iter','Algorithm','quasi-newton',...
            'MaxIterations', 1e9, 'OptimalityTolerance', eps, 'StepTolerance', 1e-6,...
            'MaxFunctionEvaluations', 1e9);
     
    else

        options                  = optimoptions('fminunc','Display','none','Algorithm','quasi-newton',...
            'MaxIterations', 1e9, 'OptimalityTolerance', eps, 'StepTolerance', 1e-6,...
            'MaxFunctionEvaluations', 1e9);
%         options                  = optimoptions('fminunc','Display','none','Algorithm','quasi-newton',...
%             'MaxIterations', 1e9, 'OptimalityTolerance', 100, ...
%             'MaxFunctionEvaluations', 1e9, 'UseParallel', Parameters.parallel);
        
    end

    %first attempt at fitting
    disp('-> Fitting a first model')

    model = vectorize_model(model, Parameters);

    f = @(x) fitwaveform_wrapper(x, model, allWfs, Parameters);
    
    [x, fout, ~] = ACD(f, length(model.vector), -0.1, 0.1, 1e4, -1e9, 1, model.vector, Parameters.solver_printout);
    %[x, ~, flag, output, del] = fminunc(f, model.vector, options);

    model.vector   = x;
    model          = devectorize_model(model, nevt, nsta, length(Parameters.t), Parameters);
    model          = evaluate(model, allWfs, Parameters, 1:nevt);

    model = vectorize_model(model, Parameters);

    disp([ '--> First solution at BIC of ' num2str(model.BIC) ]);
    rounds = 0;

    while rounds < Parameters.reset_rounds
        
        updated = false;

        modeln = model;
           
        if (rand()>0.5 && modeln.n < Parameters.max_gaussians) || model.n == 1

            modeln.n = modeln.n + 1;

            modeln.A{modeln.n} = randn()*0.1;

            a           = modeln.t0;
            a           = [a; exprnd(modeln.t(end)/4)];
            modeln.t0 = a;
            a           = modeln.w;
            a           = [a; rand()*(log(model.t(end)/20) - 2*log(model.t(2))) + 2*log(model.t(2))];
            modeln.w    = a;

            disp([' --> Lengthening RF to ' num2str(modeln.n) ]);

            modeln.vector_source = false;

        end

        if rand()>0.5 && modeln.n > 1 %shorten

            kill        = randi([ 1 modeln.n ]);
            
            a           = modeln.A{kill};
            
            if length(a) > 1
                
                kill2          = randi([ 1 length(a) ]);
                a(kill2)       = [];
                modeln.A{kill} = a;

                disp(' --> Shortening pulse order');

            else
    
                modeln.A(kill) = [];

                a           = modeln.t0;
                a(kill)     = [];
                modeln.t0   = a;
                a           = modeln.w;
                a(kill)     = [];
                modeln.w    = a;
                
                modeln.n = modeln.n - 1;

                disp([' --> Shortening RF to ' num2str(modeln.n)]);

            end

            modeln.vector_source = false;

        end

        if (rand()>0.75)

            change        = randi([ 1 modeln.n ]);

            a = modeln.A{change};
            a = [ a; randn()*0.1];
            modeln.A{change} = a;

            disp(' --> Increased a pulse order');

            modeln.vector_source = false;

        end
        
        if rand()>0.9

            disp(' --> Changing error terms');

            modeln.r = modeln.r + randn()*Parameters.r_range(2);
            modeln.f = modeln.f + randn()*Parameters.f_range(2);

            modeln.vector_source = true;

        end

        %evaluate the new model      
        modeln    = vectorize_model(modeln, Parameters);
        %search for local minima only! Hence +/- small range            
        f = @(x) fitwaveform_wrapper(x, modeln, allWfs, Parameters);
        %x = fminunc(f, modeln.vector, options);
        %[x, ~, ~] = ACD(f, length(modeln.vector), -0.1, 0.1, 500, -1e9, 1, modeln.vector, Parameters.solver_printout);

        keyboard

        modeln.vector   = x;
        modeln          = devectorize_model(modeln, nevt, nsta, length(Parameters.t), Parameters);
        modeln          = evaluate(modeln, allWfs, Parameters, 1:nevt);
    
        if (modeln.BIC + tol) < model.BIC

            model = modeln;

            model.vector   = x;
            model          = devectorize_model(model, nevt, nsta, length(Parameters.t), Parameters);
            model          = evaluate(model, allWfs, Parameters, 1:nevt);
            disp(['-> Better solution over a round found with BIC of ' num2str(model.BIC) ]);

            updated = true;

        end

        if ~updated

            rounds = rounds + 1;
            disp(['-> Round completed with no improvements. Counter raised to ' num2str(rounds) ' of ' num2str(Parameters.reset_rounds) ]);

        else

            if rounds > 0

                disp('Round completed with some improvement. Setting counter to zero');

            end
            disp('')
            rounds = 0;

        end

    end

    disp('Final RF found; exiting optimization');
    mout = model;

%     %%%%%%%%%%%%%%%%%%%%%%
%     %HMC search for errors
% 
%     if Parameters.get_errors
% 
%         f   = @(x) fitwaveform_wrapper(x, model, allWfs, Parameters, 1);
%         smp = hmcSampler(f, model.vector, 'NumSteps', 5, 'MassVectorTuningMethod', 'hessian', 'StepSize', 1e-1);
%     
%         c = [];
%     
%         disp('-> Starting HMC search for errors');
%         if Parameters.solver_printout
%     
%             smp          = tuneSampler(smp, 'VerbosityLevel', 2, 'NumPrint',1, 'NumStepSizeTuningIterations', 50, 'NumStepsLimit', 5);
%             smp.NumSteps = 5;
%     
%             for k = 1:Parameters.nchains
%     
%                 disp([ '     Chain #' num2str(k) ' of ' num2str(Parameters.nchains) ]);
%                 
%                 tmp = drawSamples(smp, 'Burnin', Parameters.burnin, 'NumSamples', Parameters.batch_size, ...
%                     'ThinSize', Parameters.batch_thin ...
%                 , 'Verbosity', 1, 'NumPrint', 1, 'StartPoint', model.vector);
%     
%                 c = [ c; tmp ];
%     
%             end
%     
%         else
%     
%             smp          = tuneSampler(smp, 'NumStepSizeTuningIterations', 50, 'NumStepsLimit', 5);
%             smp.NumSteps = 5;
%     
%             for k = 1:Parameters.nchains
%     
%                 disp([ '     Chain #' num2str(k) ' of ' num2str(Parameters.nchains) ]);
%     
%                 tmp = drawSamples(smp, 'Burnin', Parameters.burnin, 'NumSamples', Parameters.batch_size, ...
%                     'ThinSize', Parameters.batch_thin, 'StartPoint', model.vector);
%     
%                 c = [ c; tmp ];
%     
%             end
%     
%         end
%     
%         [n, ~] = size(c);
%     
%         disp('   Building models from HMC search');
%         mout = model;%save the optimimum
%         for k = 2:n+1
%     
%             m        = model;
%             m.vector = c(k-1, :)';
%             m        = devectorize_model(m, nevt, nsta, length(Parameters.t), Parameters);
%             mout(k)  = evaluate(m, allWfs, Parameters, 1:nevt, 1:nsta, []);
%     
%         end
% 
%     else
% 
%         mout = model;
% 
%     end
       
end

%     f         = @(x) fitwaveform_wrapper(x, model, allWfs, Parameters, false);
%     x = ACD( @(XV,mu) ParForParallelWrapper(f, XV, mu, 0 ), model.vector,3,1e-6,[],[],[],[],[],[],1,...
%     4,2,1,Parameters.solver_printout);

%     %search for local minima only! Hence +/- small range
%     noimprove = false;
%     best      = Inf;
% 
%     while ~noimprove
% 
%         [x, fout, ~] = ACD(f, length(model.vector), -5, 5, 1e12, -1e9, 1, model.vector, Parameters.solver_printout);
% 
%         if fout + tol < best
% 
%             best = fout;
% 
%         else
% 
%             noimprove = true;
% 
%         end
% 
%     end
%     x = ACD( @(XV,mu) ParForParallelWrapper(f, XV, mu, 0 ), model.vector,3,1e-6,[],[],[],[],[],[],1,...
%         2,4,1,Parameters.solver_printout);


%     f                               = @(x) fitwaveform_wrapper(x, model, allWfs, Parameters);
%     [x, fout, ~] = ACD(f, length(model.vector), -1, 1, 1e12, -1e9, 1, model.vector, Parameters.solver_printout);
%     x = ACD( @(XV,mu) ParForParallelWrapper(f, XV, mu, 0 ), model.vector,1,1e-6,[],[],[],[],[],[],1,...
%        2,3,2,Parameters.solver_printout);

%     exitflag = 10;
% 
%     while exitflag ~= 1
% 
%         f                               = @(x) fitwaveform_wrapper(x, model, allWfs, Parameters);
%         [x, ~, ~, fmin_out, del] = fminunc(f, model.vector, options);
% 
%         if max(abs(del)) < deltol
%             %shake it up to try to get back on track
%             %x = normrnd(x, 0.01);%nondimensionalized
%             %model.vector   = x;
%             %model          = devectorize_model(model, nevt, nsta, length(Parameters.t), Parameters, true);
% 
%             exitflag = 1;
% 
%         else
% 
%             model.vector   = x;
%             model          = devectorize_model(model, nevt, nsta, length(Parameters.t), Parameters);
% 
%         end
% 
%     end
%     f         = @(x) fitwaveform_wrapper(x, model, allWfs, Parameters, true);
%     x = ACD( @(XV,mu) ParForParallelWrapper(f, XV, mu, 0 ), model.vector,3,1e-6,[],[],[],[],[],[],1,...
%         2,2,2,Parameters.solver_printout);


%         if nsta > 1
% 
%             action = randperm(nevt + nsta + 1);%can change covarience too 
% 
%         else
% 
%             action = randperm(nevt + 1);%can change covarience too 
% 
%         end
% 
%         %do one of every transdimensional step
%         for k = 1:length(action)
% 
%             skip = true;
% 
%             modeln = model;
% 
%             index = action(k);%to avoid changing the loop index
%            
%             if index <= nevt && nsta > 1
% 
%                 addorremove = randi([ 0 1 ]);
% 
%                 if addorremove == 0 && length(modeln.sources{index}) > 2%shorten
%                     
%                     s                    = modeln.sources{index};
%                     s(end)               = [];
%                     modeln.sources{index} = s;
% 
%                     skip = false;
% 
%                     modeln.vector_source = true;
%                     disp([' --> Shortening source #' num2str(index) ' to ' num2str(length(modeln.sources{index} )) ]);
% 
%                 elseif addorremove == 1 && length(modeln.sources{index}) < length(modeln.t)
% 
%                     modeln.sources{index} = [ modeln.sources{index}; randn() ];
%                     skip = false;
% 
%                     modeln.vector_source = true;
%                     disp([' --> Lengthening source #' num2str(index) ' to ' num2str(length(modeln.sources{index} )) ]);
% 
%                 end
% 
%             end
% 
%             if nsta > 1
% 
%                 index = index - nevt;
% 
%             end
% 
%             %now loop over the station
%             for i = 1:nsta
% 
%                 addorremove = randi([ 0 1 ]);
% 
%                 if index == 1
% 
%                     if addorremove && modeln.nrf(i) > 1 %shorten
%     
%                         kill        = randi([ 1 modeln.nrf(i)]);
%                         a           = modeln.az{i};
%                         a(kill)     = [];
%                         modeln.az{i} = a;
%                         a           = modeln.Ar{i};
%                         a(kill)     = [];
%                         modeln.Ar{i} = a;
%                         a           = modeln.Br{i};
%                         a(kill)     = [];
%                         modeln.Br{i} = a;
%                         a           = modeln.At{i};
%                         a(kill)     = [];
%                         modeln.At{i} = a;
%                         a           = modeln.Bt{i};
%                         a(kill)     = [];
%                         modeln.Bt{i} = a;
%                         a           = modeln.f0{i};
%                         a(kill)     = [];
%                         modeln.f0{i} = a;
%                         a           = modeln.t0{i};
%                         a(kill)     = [];
%                         modeln.t0{i} = a;
%                         a           = modeln.w{i};
%                         a(kill)     = [];
%                         modeln.w{i} = a;
%     
%                         modeln.nrf(i) = modeln.nrf(i) - 1;
%                         disp([' --> Shortening RF at station #' num2str(i)  ' to ' num2str(modeln.nrf(i))]);
% 
%                         skip = false;
% 
%                         modeln.vector_source = false;
% 
%                     end
% 
%                     if (~addorremove && modeln.nrf(i) < Parameters.max_gaussians) || model.nrf(i) == 1
% 
%                         a           = modeln.az{i};
%                         a           = [a; randn()];
%                         modeln.az{i} = a;
%                         a           = modeln.Ar{i};
%                         a           = [a; randn()*0.1];
%                         modeln.Ar{i} = a;
%                         a           = modeln.Br{i};
%                         a           = [a; randn()*0.1];
%                         modeln.Br{i} = a;
%                         a           = modeln.At{i};
%                         a           = [a; randn()*0.1];
%                         modeln.At{i} = a;
%                         a           = modeln.Bt{i};
%                         a           = [a; randn()*0.1];
%                         modeln.Bt{i} = a;
%                         a           = modeln.t0{i};
%                         a           = [a; exprnd(modeln.t(end)/4)];
%                         modeln.t0{i} = a;
%                         a           = modeln.f0{i};
%                         a           = [a; (rand()*(maxf - minf) + minf)];
%                         modeln.f0{i} = a;
%                         a           = modeln.w{i};
%                         a           = [a; randn()];
%                         modeln.w{i} = a;
%     
%                         modeln.nrf(i) = modeln.nrf(i) + 1;
%                         disp([' --> Lengthening RF at station #' num2str(i) ' to ' num2str(modeln.nrf(i)) ]);
% 
%                         skip = false;
% 
%                         modeln.vector_source = false;
% 
%                     end
% 
%                 end
% 
%                 index = index - 1;
% 
%             end
% 
% %             if index == 1 && rand()>0.5 %still do, even if you aren't using covarience
% % 
% %                 disp(' --> Changing error terms');
% % 
% %                 modeln.r = modeln.r + randn()*Parameters.r_range(2);
% %                 modeln.f = modeln.f + randn()*Parameters.f_range(2);
% % 
% %                 skip = false;
% % 
% %                 modeln.vector_source = true;
% % 
% %             end
% 
