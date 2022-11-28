function [ model ] = blockgibbs_AB(model, allWfs, TD_parameters)
%computationally intensive!        

    %same the same form, but different realization of the random numbers.
    %Mix of wide and narrow searches. 
    pertA1 = [ exprnd(1, round(TD_parameters.line_iterations*4), 1)...
        .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1);
        exprnd(TD_parameters.step_std, round(TD_parameters.line_iterations*4), 1)...
        .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1)]*TD_parameters.ABstd;

    pertB1 = [ exprnd(1, round(TD_parameters.line_iterations*4), 1)...
        .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1);
        exprnd(TD_parameters.step_std, round(TD_parameters.line_iterations*4), 1)...
        .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1) ]*TD_parameters.ABstd;

    pertrot1 = [ exprnd(1, round(TD_parameters.line_iterations*4), 1)...
        .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1);
        exprnd(TD_parameters.step_std, round(TD_parameters.line_iterations*4), 1)...
        .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1) ]*TD_parameters.ABstd;

    %mix between small and large perturbations
    pertA1   = pertA1(randperm(length(pertA1)));
    pertB1   = pertB1(randperm(length(pertA1)));
    pertrot1 = pertrot1(randperm(length(pertrot1)));

    if TD_parameters.max_layers == 2

        pertA2 = [ exprnd(1, round(TD_parameters.line_iterations*4), 1)...
            .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1);
            exprnd(TD_parameters.step_std, round(TD_parameters.line_iterations*4), 1)...
            .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1) ]*TD_parameters.ABstd;
    
        pertB2 = [ exprnd(1, round(TD_parameters.line_iterations*4), 1)...
            .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1);
            exprnd(TD_parameters.step_std, round(TD_parameters.line_iterations*4), 1)...
            .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1) ]*TD_parameters.ABstd;

        pertrot2 = [ exprnd(1, round(TD_parameters.line_iterations*4), 1)...
            .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1);
            exprnd(TD_parameters.step_std, round(TD_parameters.line_iterations*4), 1)...
            .*(-1).^randi([0, 1], round(TD_parameters.line_iterations*4), 1) ]*TD_parameters.ABstd;

        pertA2   = pertA2(randperm(length(pertA1)));
        pertB2   = pertB2(randperm(length(pertA1)));
        pertrot2 = pertrot2(randperm(length(pertA1)));

    end

    [nevt, nsta] = size(allWfs);

    if TD_parameters.cluster

        nsearch = 1;

    else

        nsearch = nsta;

    end

    for k = 1:nsearch

        starting_lpst = model.llh/model.T + model.lp;
        
        for kk = 1:length(pertA1)

            modeln = model;
        
            modeln.A(k, 1) = model.A(k, 1) + pertA1(kk);
            modeln.B(k, 1) = model.B(k, 1) + pertB1(kk);
            modeln.fast_dir_rotation(k, 1) = model.fast_dir_rotation(k, 1) + pertrot1(kk);

            if TD_parameters.max_layers == 2
    
                modeln.A(k, 2) = model.A(k, 2) + pertA2(kk);
                modeln.B(k, 2) = model.B(k, 2) + pertB2(kk);
                modeln.fast_dir_rotation(k, 2) = model.fast_dir_rotation(k, 2) + pertrot2(kk);

            end

            modeln = make_dtphi(modeln);

            if TD_parameters.cluster

                modeln = evaluate(modeln, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
        
            else

                modeln = evaluate(modeln, allWfs, TD_parameters, 1:nevt, k, []);

            end

            lpst(kk) = modeln.llh/model.T + modeln.lp;

        end

        pA1 = [ pertA1;   0];
        pB1 = [ pertB1;   0];
        pr1 = [ pertrot1; 0];
        pA2 = [ pertA2;   0];
        pB2 = [ pertB2;   0];
        pr2 = [ pertrot2; 0];

        if ~TD_parameters.optimize
    
            %calculate relative probablities and make a random selection

            pst  = [ lpst starting_lpst];
            pst  = exp(pst-max(pst));%avoid overflow. Underflow ok (just rounds to zero)
            pst  = pst/sum(pst);%make it a probability
            i    = randsample(1:length(pA2), 1, 'true', pst);
    
            model.A(k, 1)                 = model.A(k, 1) + pA1(i);
            model.B(k, 1)                 = model.B(k, 1) + pB1(i);
            model.fast_dir_rotation(k, 1) = model.fast_dir_rotation(k, 1) + pr1(i);

            if TD_parameters.max_layers == 2
    
                model.A(k, 2) = model.A(k, 2) + pA2(i);
                model.B(k, 2) = model.B(k, 2) + pB2(i);
                model.fast_dir_rotation(k, 2) = model.fast_dir_rotation(k, 2) + pr2(i);

            end
    
        else
    
            [mpst, i] = max(lpst);
    
            if mpst > starting_lpst
    
                model.A(k, 1) = model.A(k, 1) + pertA1(i);
                model.B(k, 1) = model.B(k, 1) + pertB1(i);
                model.fast_dir_rotation(k, 1) = model.fast_dir_rotation(k, 1) + pr1(i);
    
                if TD_parameters.max_layers == 2
        
                    model.A(k, 2) = model.A(k, 2) + pertA2(i);
                    model.B(k, 2) = model.B(k, 2) + pertB2(i);
                    model.fast_dir_rotation(k, 2) = model.fast_dir_rotation(k, 2) + pr2(i);

                end
        
            end
    
        end

        model        = make_dtphi(model);

        if TD_parameters.cluster
    
            model = evaluate(model, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
    
        else
    
            model = evaluate(model, allWfs, TD_parameters, 1:nevt, k, []);
    
        end

    end

    %model = vectorize_model(model, TD_parameters);
    
end