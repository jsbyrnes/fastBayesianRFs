function [BIC, grad] = fitwaveform_wrapper(vector, model, allWfs, Parameters)

    [nevt, nsta] = size(allWfs);

    model.vector = vector;
    model        = devectorize_model(model, nevt, nsta, length(Parameters.t), Parameters);
    
    model = evaluate(model, allWfs, Parameters, 1:nevt);
    BIC   = model.BIC;%fminunc is minimization

    if nargout > 1 % gradient required

        model = fill_gradient_waveforms(model, allWfs, Parameters);
        grad  = -2*model.del';%AIC is the objective

    end

end
