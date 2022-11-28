function [rf_mean, rf_std, mout ] = fit_waveforms_rjMCMC(allWfs, Parameters)

    mout = [];

    sig = Parameters.sig;

    total         = 0;
    acceptedtotal = 0;

    %minf    = 1/(Parameters.t(end));
    %maxf    = 0.5/(Parameters.t(2) - Parameters.t(1));

    [nevt, nsta] = size(allWfs);

    disp([ '-> Building a starting models for ' allWfs.station ])

    for k = 1:Parameters.nchains

        model_set(k) = make_starting_models(Parameters, allWfs);

        if k <= ceil(Parameters.nchains*0.2)

            model_set(k).temperature = 1;

        end

    end
    %model = vectorize_model(model, Parameters);

    for iter = 1:Parameters.niterations

        if mod(iter, Parameters.printon) == 0

            disp([ allWfs.station ' at #' num2str(iter) ' with lpst of ' num2str(mean([model_set([model_set(:).temperature]==1).lpst])) ]);

        end

        if mod(iter, Parameters.update_sig) && iter < Parameters.burnin

            if acceptedtotal/total > 0.2

                sig = sig/0.75;

            elseif acceptedtotal/total < 0.1

                sig = sig*0.75;
                
            end

            acceptedtotal = 0;
            total         = 0;

        end

        for k = 1:Parameters.nchains

            model = model_set(k);
            
            accepted = false;
            action   = randi([1 8]);
        
            switch action 
    
                case 1 %shape/height
    
                    modeln = model;

                    if model.n > 1

                        index  = randi([ 1 model.n]);

                    else

                        index = 1;

                    end
                        
                    A = modeln.A{index};
    
                    A = A + 0.1*randn(size(A))*sig;
    
                    modeln.A{index} = A;
    
                    modeln.update_rf = true;

                    modeln = evaluate(modeln, allWfs, Parameters);
    
                    if log(rand()) < ((modeln.lpst - model.lpst)/model.temperature)
    
                        model = modeln;
                        accepted = true;

                    end
    
                case 2 %time
    
                    modeln = model;
                    if model.n > 1

                        index  = randi([ 1 model.n]);

                    else

                        index = 1;

                    end
                            
                    modeln.t0(index) = modeln.t0(index) + 0.2*randn()*sig*model.t(end);%time is a smaller step
        
                    modeln.update_rf = true;

                    modeln = evaluate(modeln, allWfs, Parameters);
    
                    if log(rand()) < ((modeln.lpst - model.lpst)/model.temperature)
    
                        model = modeln;
                        accepted = true;

                    end
    
                case 3 %width
    
                    modeln = model;
                    if model.n > 1

                        index  = randi([ 1 model.n]);

                    else

                        index = 1;

                    end
                            
                    modeln.w(index)  = modeln.w(index) + randn()*sig;
        
                    modeln.update_rf = true;

                    modeln = evaluate(modeln, allWfs, Parameters);
    
                    if log(rand()) < ((modeln.lpst - model.lpst)/model.temperature)
    
                        model = modeln;
                        accepted = true;

                    end
     
                case 4 %cov
    
                    modeln = model;
                        
                    modeln.r = modeln.r + randn()*Parameters.r_range(2)*sig;
                    modeln.f = modeln.f + randn()*Parameters.f_range(2)*sig;
                    modeln   = build_C(modeln);
                    modeln.update_rf = false;
                    modeln   = evaluate(modeln, allWfs, Parameters);
    
                    if log(rand()) < ((modeln.lpst - model.lpst)/model.temperature)
    
                        model = modeln;
                        accepted = true;

                    end
    
                case 5 %birth
    
                    modeln = model;
    
                    birth = randi([ 1 (modeln.n+1) ]);

                    if birth == (modeln.n+1)

                        modeln.n = modeln.n + 1;
            
                        modeln.A{modeln.n} = randn()*0.1;
            
                        a           = modeln.t0;
                        a           = [a; exprnd(modeln.t(end)/4)];
                        modeln.t0 = a;
                        a           = modeln.w;
                        a           = [a; rand()*(log(model.t(end)/20) - 2*log(model.t(2))) + 2*log(model.t(2))];
                        modeln.w    = a;

                    else

                        a = modeln.A{birth};
                        a = [ a; randn()*0.1];
                        modeln.A{birth} = a;

                    end
    
                    modeln.update_rf = true;
                    modeln   = evaluate(modeln, allWfs, Parameters);
    
                    if log(rand()) < ((modeln.lpst - model.lpst)/model.temperature)
    
                        model = modeln;
                        accepted = true;

                    end
    
                case 6 %death
    
                    modeln = model;

                    if model.n > 1

                        kill  = randi([ 1 model.n]);

                    else

                        kill = 1;

                    end
            
                    a           = modeln.A{kill};
                    
                    if length(a) > 1
                        
                        kill2          = randi([ 1 length(a) ]);
                        a(kill2)       = [];
                        modeln.A{kill} = a;
                
                    elseif length(a)==1 && model.n > 1
            
                        modeln.A(kill) = [];
        
                        a           = modeln.t0;
                        a(kill)     = [];
                        modeln.t0   = a;
                        a           = modeln.w;
                        a(kill)     = [];
                        modeln.w    = a;
                        
                        modeln.n = modeln.n - 1;
                
                    end

                    modeln.update_rf = true;

                    modeln   = evaluate(modeln, allWfs, Parameters);
    
                    if log(rand()) < ((modeln.lpst - model.lpst)/model.temperature)
    
                        model = modeln;
                        accepted = true;

                    end
    
                case 7 %sig pairs
    
                    modeln = model;
        
                    modeln.sig = modeln.sig + randn(size(modeln.sig))*Parameters.sig_range(2)*sig;
                    modeln.update_rf = false;
                    modeln   = evaluate(modeln, allWfs, Parameters);
    
                    if log(rand()) < ((modeln.lpst - model.lpst)/model.temperature)
    
                        model = modeln;
                        accepted = true;

                    end
    
            end

            if model.temperature == 1 && action ~= 5 && action ~= 6

                total = total + 1;

                if accepted

                    acceptedtotal = acceptedtotal + 1;

                end

            end

            model_set(k) = model;

        end

        inds = randperm(Parameters.nchains);

        for k1 = 1:length(inds)

            for k2 = 1:length(inds)

                if k2>=k1

                    continue

                end

                alpha = ((model_set(k2).lpst - model_set(k1).lpst)/model_set(k1).temperature ...
                                + (model_set(k1).lpst - model_set(k2).lpst)/model_set(k2).temperature);

                if alpha > log(rand())

                    T1 = model_set(k1).temperature;
                    model_set(k1).temperature = model_set(k2).temperature;
                    model_set(k2).temperature = T1;

                end

            end
            
        end

        if mod(iter, Parameters.saveint)==0

            if iter > Parameters.burnin

                mout = [ mout model_set([model_set(:).temperature]==1) ];

            end

        end

        history(iter) = mean([model_set([model_set(:).temperature]==1).lpst]);

    end

    %build the RF
    rf_mean = zeros(size(model.t));
    rf_std  = zeros(size(model.t));

    for k = 1:length(mout)

        rf_matrix(k, :) = mout(k).rf;

    end

    rf_mean = mean(rf_matrix);
    rf_std  = std(rf_matrix);

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
