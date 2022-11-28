function model = neldermead(model, allWfs, TD_parameters)

    alpha = 1;
    gamma = 2;
    rho   = 0.5;
    sigma = 0.5;

    [nevt, nsta] = size(allWfs);
    n            = length(TD_parameters.t);

    disp('Making the simplex')
    model = vectorize_model(model, TD_parameters);
    %simplex(1) = model;
    for k = 1:(model.nparam + 1)

        simplex(k) = make_starting_models_NM(model, TD_parameters, allWfs);

    end

    maxpert = Inf;

    niter = 0;

    [~, ind] = sort([simplex(:).lpst], 'descend');
    simplex  = simplex(ind);

    while rms(simplex(end).vector - simplex(1).vector) > 1e-4

        [~, ind] = sort([simplex(:).lpst], 'descend');
        simplex  = simplex(ind);

        niter = niter + 1;

        if mod(niter, 10) == 1

            disp([ 'Iteration ' num2str(niter) ' at peak llh ' num2str(simplex(1).llh, 5) ...
                ' mean llh of ' num2str(mean([simplex(:).llh]), 5) ]);
        end

        hist_llh(niter) = simplex(1).llh;
        hist_lp(niter)  = simplex(1).lp;

        xmean = centroid(simplex);
        xr    = xmean + alpha*(xmean - simplex(end).vector);

        %reflection
        mr        = simplex(end);
        mr.vector = xr;
        mr        = devectorize_model(mr, nevt, nsta, n, TD_parameters);
        mr        = evaluate(mr, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

        if mr.lpst < simplex(1).lpst && mr.lpst > simplex(end-1).lpst

            simplex(end) = mr;
            %disp('Reflection')
            continue
        
        else

            if mr.lpst > simplex(1).lpst

                %expansion
                xe        = mr.vector + gamma*(mr.vector - xmean);
                me        = simplex(end);
                me.vector = xe;
                me        = devectorize_model(me, nevt, nsta, n, TD_parameters);
                me        = evaluate(me, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

                if me.lpst > mr.lpst

                    simplex(end) = me;

                else

                    simplex(end) = mr;

                end

                disp('Expansion')
                continue

            else

                %contraction
                if mr.lpst >= simplex(end).lpst

                    xc        = xmean + rho*(mr.vector - xmean);
                    mc        = simplex(end);
                    mc.vector = xc;
                    mc        = devectorize_model(mc, nevt, nsta, n, TD_parameters);
                    mc        = evaluate(mc, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

                    if mc.lpst > mr.lpst

                        simplex(end) = mc;
                        %disp('Contraction')
                        continue

                    end

                else

                    xc        = xmean + rho*(simplex(end).vector - xmean);
                    mc        = simplex(end);
                    mc.vector = xc;
                    mc        = devectorize_model(mc, nevt, nsta, n, TD_parameters);
                    mc        = evaluate(mc, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

                    if mc.lpst > simplex(end).lpst

                        simplex(end) = mc;
                        %disp('Contraction')
                        continue

                    end

                end

                %shrinkage
                disp('Shrinking')
                for k = 2:(model.nparam + 1)
            
                    xs         = simplex(1).vector + sigma*(simplex(k).vector - simplex(1).vector);
                    ms         = simplex(k);
                    ms.vector  = xs;
                    ms         = devectorize_model(ms, nevt, nsta, n, TD_parameters);
                    simplex(k) = evaluate(ms, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
            
                end

            end

        end

    end

    [~, ind] = sort([simplex(:).lpst], 'descend');
    simplex  = simplex(ind);
    model    = simplex(1);

end

function xmean = centroid(simplex)
    %assume

    [~, ind] = sort([simplex(:).lpst], 'descend');
    simplex  = simplex(ind);

    n = length(simplex)-1;

    xmean = simplex(1).vector;

    for k = 2:n
        xmean = xmean + simplex(k).vector;
    end

    xmean = xmean/n;

end