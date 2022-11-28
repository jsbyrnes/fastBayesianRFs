function model = devectorize_model(model, nevt, nsta, n, Parameters)

    %you act on the model vector, this fills it back out for evaluation
    %source parameters

    vec = model.vector;%will clip for simplicity
    dt  = model.t(2) - model.t(1);

    for k = 1:model.n

        model.A{k} = vec(1:length(model.A{k}))*0.1;
        vec(1:length(model.A{k})) = [];

    end

    model.t0 = vec(1:model.n);
    vec(1:model.n) = [];
    model.w  = vec(1:model.n);
    vec(1:model.n) = [];

    %error terms
    model.sig = reshape((vec(1:nevt*nsta)*Parameters.sig_range(2)) + Parameters.sig_range(1), [ nevt nsta ]);
    vec = vec((nevt*nsta + 1):end);
    %model.sigT = reshape((vec(1:nevt*nsta)*Parameters.sig_range(2)) + Parameters.sig_range(1), [ nevt nsta ]);
    %vec = vec((nevt*nsta + 1):end);

    if Parameters.use_covarience

        model.r = vec(1)*Parameters.r_range(2) + Parameters.r_range(1);
        vec(1)  = [];
        model.f = vec(1:length(model.f))*Parameters.f_range(2) + Parameters.f_range(1);
        vec(1:length(model.f)) = [];
        %model.a = vec(1:length(model.a));
        %vec(1:length(model.a)) = [];

        %now apply parameters to useful model
        model = build_C(model);

    else

        model.Cinv   = eye(n);
        model.logdet = 1;
        
    end

end