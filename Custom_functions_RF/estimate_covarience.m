function model = estimate_covarience(TD_parameters, allWfs, model)

    c = zeros(length(TD_parameters.t), 1);

    d = [];

    for k = 1:numel(allWfs)

        m = max( abs(allWfs(k).north));
        s = std( allWfs(k).north);

        if kurtosis(model.N{k} - allWfs(k).north) < 8
        %if m < TD_parameters.noisy_threshold*s

            ctmp = ifft(fft(model.N{k} - allWfs(k).north).*conj(fft(model.N{k} - allWfs(k).north)));
            ctmp = ctmp/ctmp(1);

            c = c + ctmp;

            d = [ d allWfs(k).north/rms(allWfs(k).north) ];

        end

        m = max( abs(allWfs(k).east));
        s = std( allWfs(k).east);

        if kurtosis(model.E{k} - allWfs(k).east) < 8
        %if m < TD_parameters.noisy_threshold*s

            ctmp = ifft(fft(model.E{k} - allWfs(k).east).*conj(fft(model.E{k} - allWfs(k).east)));
            ctmp = ctmp/ctmp(1);

            c = c + ctmp;

            d = [ d allWfs(k).east/rms(allWfs(k).east) ];

        end

    end

    disp([ num2str(c(1)) ' traces used for initial covarience estimate' ])

    c = c/c(1);

    c = c(1:round(length(TD_parameters.t)/2));

    %now find a best fit

    dt = TD_parameters.t(1:round(length(TD_parameters.t)/2)) - TD_parameters.t(1);
    fit = @(x) sum((c - exp(-exp(x(1))*dt).*cos(2*pi*exp(x(2))*dt)).^2);

    options                  = optimoptions('fminunc','Display','none','Algorithm','quasi-newton',...
        'MaxIterations', 1e9, 'OptimalityTolerance', 1e-4, 'FiniteDifferenceType', 'central', ...
        'MaxFunctionEvaluations', 1e9);

    [x, feval, exitflag, fmin_out, del] = fminunc( fit, [model.r model.f], options);
    %x = fminunc( fit, [model.r model.f], options);

    model.r = x(1);
    model.f = x(2);

    model = build_C(model, TD_parameters.t);

end


%     %first, do a grid search
%     [R, F] = meshgrid( linspace(TD_parameters.r_range(1) - TD_parameters.r_range(2)*4, ...
%         TD_parameters.r_range(1) + TD_parameters.r_range(2)*4, 100), ...
%         linspace(TD_parameters.f_range(1) - TD_parameters.f_range(2)*4, ...
%         log(TD_parameters.high_pass*2), 100));
% 
%     for k = 1:numel(R)
% 
%         fit_search(k) = fit( [ R(k), F(k) ]);
% 
%     end

    %contourf(R, F, reshape(log10(fit_search), size(R)), -1:0.05:0)

    %then refine
    %[~, ind] = min(fit_search);

