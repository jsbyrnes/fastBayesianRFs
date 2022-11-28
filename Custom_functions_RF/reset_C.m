function [ model ] = reset_C(model, allWfs, TD_parameters)
        
    %number of samples scale to however many times you update these 

    perts1 = [ exprnd(TD_parameters.r_range(2), round(TD_parameters.line_iterations), 1)...
            .*(-1).^randi([0, 1], round(TD_parameters.line_iterations), 1);
            exprnd(TD_parameters.step_std*TD_parameters.r_range(2), round(TD_parameters.line_iterations), 1)...
            .*(-1).^randi([0, 1], round(TD_parameters.line_iterations), 1)];

    perts2 = [ exprnd(TD_parameters.f_range(2), round(TD_parameters.line_iterations), 1)...
            .*(-1).^randi([0, 1], round(TD_parameters.line_iterations), 1);
            exprnd(TD_parameters.step_std*TD_parameters.f_range(2), round(TD_parameters.line_iterations), 1)...
            .*(-1).^randi([0, 1], round(TD_parameters.line_iterations), 1)];

    [nevt, nsta] = size(allWfs);

    for kk = 1:length(perts1)

        modeln = model;

        modeln.r         = model.r + perts1(kk);
        modeln.f         = model.f + perts2(kk);
        modeln           = build_C(modeln, TD_parameters.t);
        modeln           = evaluate(modeln, allWfs, TD_parameters, 1:nevt, 1:nsta, 'error');

        lps(kk) = modeln.llh + modeln.lp;

    end

    [ml, i] = max(lps);

    if ml > (model.llh + model.lp)

        modeln.r         = model.r + perts1(i);
        modeln.f         = model.f + perts2(i);
        modeln           = build_C(modeln, TD_parameters.t);
        model            = evaluate(modeln, allWfs, TD_parameters, 1:nevt, 1:nsta, 'error');%burn a run

    end

    model = vectorize_model(model, TD_parameters);
    
end