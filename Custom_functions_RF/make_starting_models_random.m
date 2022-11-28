function model = make_starting_models_random(TD_parameters, allWfs)

    t = TD_parameters.t;

    [nevt, nsta] = size(allWfs);

    model.nparam = length(TD_parameters.t)*nevt + nevt + nevt*nsta*4 + nsta*3*TD_parameters.max_layers;

    if ~TD_parameters.orientations

        model.sta_or            = zeros(nsta,1)*TD_parameters.orientation_std;%approximation to von mises

    else

        model.nparam = model.nparam + nsta;

    end

    if ~TD_parameters.covariences 
        
        model.r = TD_parameters.r_range(1);% + randn()*TD_parameters.r_range(2);%rand(n,1)*TD_parameters.max_dt;%log(0.9);%
        model.f = TD_parameters.f_range(1);% + randn()*TD_parameters.f_range(2);%rand(n,1)*TD_parameters.max_dt;%log(0.5);

    else

        model.nparam = model.nparam + 2;

    end

    model.vector        = randn(model.nparam, 1);
    model               = devectorize_model(model, nevt, nsta, length(TD_parameters.t), TD_parameters);
    model.polarization0 = TD_parameters.polarization;

    model = build_C(model, t);
    model = evaluate(model, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

    model.id = randi([ 1 1e4 ]);

end

% n = floor(length(TD_parameters.t)/3);
% 
% for k = 1:(nevt*nsta)
% 
%     %model.amp(k) = log(max([ allWfs(k).north; allWfs(k).east ]));%peak is real peak, scaled down to damp gradient
%     %model.sig(k) = log(std([ allWfs(k).north(1:n); allWfs(k).east(1:n) ])/4);%~quarter variation is noise
% 
% end
% 
% model.amp = reshape(model.amp, size(allWfs));
% %model.sig = reshape(model.sig, size(allWfs));
