function model = make_starting_models_GA(Parameters, allWfs)

    t = Parameters.t;

    [nevt, nsta] = size(allWfs);

    if nsta > 1 && Parameters.wavelet

        %per evt parameters
        for k = 1:nevt
        
           wavelet = zeros(length(Parameters.t),1);
    
            %stack the stack to the nth root
            for j = 1:nsta
    
                %rotate to radial
                radial = allWfs(k,j).north*cos(Parameters.polarization(k)) + ...
                    allWfs(k,j).east*sin(Parameters.polarization(k));
    
                if ~any(isnan(radial))
    
                    wavelet = wavelet + sign(radial).*abs(radial).^(1/Parameters.nth_root);
    
                end
    
            end
    
            model.wavelet(k, :) = sign(wavelet).*abs(wavelet).^(Parameters.nth_root);
            model.wavelet(k, :) = (model.wavelet(k, :)/rms(model.wavelet(k, :)) + 0.25*randn(size(model.wavelet(k, :))));
            %model.wavelet(k, :) = randn(length(Parameters.t),1);%model.wavelet(k, :)/rms(model.wavelet(k, :));% + randn();
    
            if Parameters.use_polarization
    
                model.polarization(k,1) = Parameters.polarization_std*randn();
    
            else
    
                model.polarization(k,1) = 0;
    
            end
    
        end

    end

    model.polarization0 = Parameters.polarization;

   for k = 1:nsta

        x = [allWfs(:, k).orientation_k];

        if isempty(x(~isnan(x)))
            sta_err(k,1) = 2*pi;
        else
            sta_err(k,1) = unique(x(~isnan(x)));
            sta_err(k,1) = sqrt(1/sta_err(k,1));
        end

    end

    if Parameters.prior_information.use

        %per station parameters
        if Parameters.cluster
    
            model.dt(1, 1) = Parameters.prior_information.dt(1) + Parameters.prior_information.dt_std(1)*randn();
            model.dt(1, 2) = Parameters.prior_information.dt(2) + Parameters.prior_information.dt_std(2)*randn();
            model.fast_dir(1,1) = Parameters.prior_information.phi(1) + Parameters.prior_information.phi_std(1)*randn();%in radians
            model.fast_dir(1,2) = Parameters.prior_information.phi(2) + Parameters.prior_information.phi_std(2)*randn();%in radians
            model   = make_AB(model);
            model   = make_dtphi(model);%enforce a pi/2 wrap
            model.fast_dir_rotation(1,1) = Parameters.prior_information.rot(1) + Parameters.prior_information.rot_std(2)*randn();%in radians
            model.fast_dir_rotation(1,2) = Parameters.prior_information.rot(2) + Parameters.prior_information.rot_std(2)*randn();%in radians

%              model.dt(1, 1) = Parameters.prior_information.dt(1);
%              model.dt(1, 2) = Parameters.prior_information.dt(2);
%              model.fast_dir(1,1) = Parameters.prior_information.phi(1);%in radians
%              model.fast_dir(1,2) = Parameters.prior_information.phi(2);%in radians
%              model   = make_AB(model);
%              model   = make_dtphi(model);%enforce a pi/2 wrap
%              model.fast_dir_rotation(1,1) = Parameters.prior_information.rot(1);%in radians
%              model.fast_dir_rotation(1,2) = Parameters.prior_information.rot(2);%in radians

        else
    
            model.dt(:, 1) = randn(nsta,1)*Parameters.prior_information.dt_std(1) + Parameters.prior_information.dt(1);
            model.dt(:, 2) = randn(nsta,1)*Parameters.prior_information.dt_std(2) + Parameters.prior_information.dt(2);
            model.fast_dir(:,1) = randn(nsta,1)*Parameters.prior_information.phi_std(1) + Parameters.prior_information.phi(1);%in radians
            model.fast_dir(:,2) = randn(nsta,1)*Parameters.prior_information.phi_std(2) + Parameters.prior_information.phi(2);%in radians
            model   = make_AB(model);
            model   = make_dtphi(model);%enforce a pi/2 wrap
            model.fast_dir_rotation(:,1) = randn(nsta,1)*Parameters.prior_information.rot_std(1) + Parameters.prior_information.rot(1);%in radians
            model.fast_dir_rotation(:,2) = randn(nsta,1)*Parameters.prior_information.rot_std(2) + Parameters.prior_information.rot(2);%in radians
    
        end

    else

        %per station parameters

        %note the sqrt(2) since dt = sqrt(A^2 + B^2)
        
        if Parameters.cluster
    
            model.A = randn(1,Parameters.max_layers)*Parameters.dtstd/sqrt(2);
            model.B = randn(1,Parameters.max_layers)*Parameters.dtstd/sqrt(2);
            model   = make_dtphi(model);
            model.fast_dir_rotation = randn(1,Parameters.max_layers)*Parameters.rotation_std;%in radians
        
        else
    
            model.A = randn(nsta,Parameters.max_layers)*Parameters.dtstd/sqrt(2);
            model.B = randn(nsta,Parameters.max_layers)*Parameters.dtstd/sqrt(2);
            model   = make_dtphi(model);
            model.fast_dir_rotation = randn(nsta,Parameters.max_layers)*Parameters.rotation_std;%in radians
    
        end

    end

    if Parameters.use_tSA

        if Parameters.cluster
    
            model.tSA = randn(1,Parameters.max_layers)*Parameters.tSAstd;
        
        else
    
            model.tSA = randn(nsta,Parameters.max_layers)*Parameters.tSAstd;
    
        end

    else

        if Parameters.cluster
    
            model.tSA = zeros(1,Parameters.max_layers);
        
        else
    
            model.tSA = zeros(nsta,Parameters.max_layers);
    
        end
        
    end

    if Parameters.use_orientations

        model.sta_or = sta_err*randn();

    else

        model.sta_or = zeros(nsta,1);%approximation to von mises

    end

    %per station-evt pairs
    model.sig      = zeros(nevt, nsta);

    if nsta > 1 && Parameters.wavelet

        model.amp      = zeros(nevt, nsta);%log
        model.shift    = zeros(nevt, nsta);
        model.dtS      = Parameters.dtS(1)*ones(nevt, nsta);%*Parameters.max_dtS;%ones(20, 1);% abs(randn(n,1)*Parameters.dt_std);

    end

    %global parameters

    if Parameters.use_covarience

        model.r = Parameters.r_range(1) + randn()*Parameters.r_range(2);%rand(n,1)*Parameters.max_dt;%log(0.9);%
        model.f = Parameters.f_range(1) + randn()*Parameters.f_range(2);%rand(n,1)*Parameters.max_dt;%log(0.5);

    else
        
        model.r = 0;
        model.f = 0;

    end

    model.N = cellfun(@ (x) zeros(length(Parameters.t),1), cell(nevt, nsta), 'UniformOutput', false);
    model.E = cellfun(@ (x) zeros(length(Parameters.t),1), cell(nevt, nsta), 'UniformOutput', false);
    model   = build_C(model, t);
    model   = evaluate(model, allWfs, Parameters, 1:nevt, 1:nsta, []);

    %vectorize the model
    model             = vectorize_model(model, Parameters);
    
    model.id    = randi([ 1 1e4 ]);

end
