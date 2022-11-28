function [ model, iter, pert ] = newton_C(model, allWfs, TD_parameters)
        
    [nevt, nsta] = size(allWfs);

    h = 1e-3;
    hr = h;%TD_parameters.r_range(2)*TD_parameters.h;
    hf = h;%TD_parameters.f_range(2)*TD_parameters.h;

    step = 1e9;
    tol  = h;

    iter = 0;

    pert = 1e9;

    while (pert > tol)

        r0 = model.r;
        f0 = model.f;

        %first derivatives
        modelpr   = model;
        modelpr.r = model.r + hr;
        modelpr   = apply_update(modelpr, 'r', TD_parameters.t);
        modelpr   = evaluate(modelpr, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
    
        modelnr   = model;
        modelnr.r = model.r - hr;
        modelnr   = apply_update(modelnr, 'r', TD_parameters.t);
        modelnr   = evaluate(modelnr, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
    
        modelpf   = model;
        modelpf.f = model.f + hf;
        modelpf   = apply_update(modelpf, 'r', TD_parameters.t);
        modelpf   = evaluate(modelpf, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
    
        modelnf   = model;
        modelnf.f = model.f - hf;
        modelnf   = apply_update(modelnf, 'r', TD_parameters.t);
        modelnf   = evaluate(modelnf, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
    
        %second derivatives
        modelprpf   = model;
        modelprpf.r = model.r + hr;
        modelprpf.f = model.f + hf;
        modelprpf   = apply_update(modelprpf, 'r', TD_parameters.t);
        modelprpf   = evaluate(modelprpf, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
    
        modelprnf   = model;
        modelprnf.r = model.r + hr;
        modelprnf.f = model.f - hf;
        modelprnf   = apply_update(modelprnf, 'r', TD_parameters.t);
        modelprnf   = evaluate(modelprnf, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
    
        modelnrpf   = model;
        modelnrpf.r = model.r - hr;
        modelnrpf.f = model.f + hf;
        modelnrpf   = apply_update(modelnrpf, 'r', TD_parameters.t);
        modelnrpf   = evaluate(modelnrpf, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
    
        modelnrnf   = model;
        modelnrnf.r = model.r - hr;
        modelnrnf.f = model.f - hf;
        modelnrnf   = apply_update(modelnrnf, 'r', TD_parameters.t);
        modelnrnf   = evaluate(modelnrnf, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

        d1(1,1) = (modelpr.llh + modelpr.lp - modelnr.llh - modelnr.lp)/(2*hr);
        d1(2,1) = (modelpf.llh + modelpf.lp - modelnf.llh - modelnf.lp)/(2*hf);
    
        H(1,1) = ((modelpr.llh + modelpr.lp) + (modelnr.llh + modelnr.lp) - 2*(model.llh + model.lp))/hr^2;
        H(2,2) = ((modelpf.llh + modelpf.lp) + (modelnf.llh + modelnf.lp) - 2*(model.llh + model.lp))/hf^2;
        H(1,2) = ( (modelprpf.llh + modelprpf.lp) + (modelnrnf.llh + modelnrnf.lp) ...
            - (modelprnf.llh + modelprnf.lp) - (modelnrpf.llh + modelnrpf.lp))/(4*hr*hf);
        H(2,1) = H(1,2);
    
        step = H\d1;
    
%         bunch = -1:0.05:0.25;
%         for k = 1:length(bunch)
%             modelb   = model;
%             modelb.r = model.r - bunch(k);
%             %modelb.f = model.f - hf;
%             modelb   = apply_update(modelb, 'r', TD_parameters.t);
%             modelb   = evaluate(modelb, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
% 
%             results(k) = modelb.lpst;
% 
%         end
%         plot(bunch, results)

        modeln   = model;
        modeln.r = modeln.r - step(1);
        modeln.f = modeln.f - step(2);
        modeln   = apply_update(modeln, 'r', TD_parameters.t);
        modeln   = evaluate(modeln, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
    
        if (modeln.llh + modeln.lp) > (model.llh + model.lp) && sum(diag(H)<0)==2
    
            model = modeln;
    
        else
    
            alpha   = 1e-3;
            step(:) = 1e9;

            niter   = 0;

            %may terminate at higher misfit, trying to advance anyway.
            %Probably means you have positive curvatures
            while (modeln.lpst <= model.lpst) && sum(abs(step)>(5*h))==2
    
                step(1)  = alpha*d1(1);
                step(2)  = alpha*d1(2);
                modeln   = model;
                modeln.r = modeln.r + step(1);
                modeln.f = modeln.f + step(2);
                modeln   = apply_update(modeln, 'r', TD_parameters.t);
                modeln   = evaluate(modeln, allWfs, TD_parameters, 1:nevt, 1:nsta, []);
    
                alpha = alpha/2;                        
                niter = niter + 1;
                
            end
                    
            model = modeln;
        
        end

        pert = max([ abs(model.r - r0) abs(model.f - f0) ]);
        iter = iter + 1;

    end

end