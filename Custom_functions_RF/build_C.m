function model = build_C(model)

    t = model.t;

    C = zeros(length(t), length(t));
    
    f = exp(model.f);
    r = exp(model.r);

    for i = 1:length(t)
    
        for j = 1:length(t)
    
            dt     = abs(t(j) - t(i));
            C(i,j) = exp(-r*dt);
    
            for k = 1:length(f)

                C(i,j) = C(i,j)*cos(2*pi*f(k)*dt);

            end

        end
    
    end

    try

        model.Cinv = inv(C);
        model.logdet   = logdet(C, 'chol');

    catch

        model.Cinv   = eye(length(t));
        model.logdet = 1e9;%occasionally occurs, r is too large. Reject the model 
    
    end

end