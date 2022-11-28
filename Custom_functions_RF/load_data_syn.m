function [ allWfs, synmodel, TD_parameters ] = load_data_syn(TD_parameters)

    nevt = 5;
    nsta = 10;

    synmodel.wavelet(1, :) = exp(-((TD_parameters.t - mean(TD_parameters.t)).^2)/(2*3.^2));
    synmodel.wavelet(1, :) = [ 0 diff(synmodel.wavelet(1, :))/(TD_parameters.t(2) - TD_parameters.t(1)) ];
    synmodel.wavelet(1, :) = (synmodel.wavelet(1, :)/rms(synmodel.wavelet(1, :)));

    synmodel.wavelet(2, :) = -exp(-((TD_parameters.t - mean(TD_parameters.t)).^2)/(2*3.^2));
    synmodel.wavelet(2, :) = [ 0 diff(synmodel.wavelet(2, :))/(TD_parameters.t(2) - TD_parameters.t(1)) ];
    synmodel.wavelet(2, :) = [ 0 diff(synmodel.wavelet(2, :))/(TD_parameters.t(2) - TD_parameters.t(1)) ];%second derivative
    synmodel.wavelet(2, :) = (synmodel.wavelet(2, :)/rms(synmodel.wavelet(2, :)));

    synmodel.wavelet(3, :) = exp(-((TD_parameters.t - mean(TD_parameters.t) - 6).^2)/(2*3.^2)) ...
        + exp(-((TD_parameters.t - mean(TD_parameters.t) + 6).^2)/(2*3.^2));
    synmodel.wavelet(3, :) = [ 0 diff(synmodel.wavelet(3, :))/(TD_parameters.t(2) - TD_parameters.t(1)) ];
    synmodel.wavelet(3, :) = (synmodel.wavelet(3, :)/rms(synmodel.wavelet(3, :)));

    synmodel.wavelet(4, :) = exp(-((TD_parameters.t - mean(TD_parameters.t) - 8).^2)/(2*5.^2)) ...
        + exp(-((TD_parameters.t - mean(TD_parameters.t) + 8).^2)/(2*2.^2));
    synmodel.wavelet(4, :) = [ 0 diff(synmodel.wavelet(4, :))/(TD_parameters.t(2) - TD_parameters.t(1)) ];
    synmodel.wavelet(4, :) = (synmodel.wavelet(4, :)/rms(synmodel.wavelet(4, :)));

    synmodel.wavelet(5, :) = -exp(-((TD_parameters.t - mean(TD_parameters.t)).^2)/(2*8.^2));
    synmodel.wavelet(5, :) = [ 0 diff(synmodel.wavelet(5, :))/(TD_parameters.t(2) - TD_parameters.t(1)) ];
    synmodel.wavelet(5, :) = [ 0 diff(synmodel.wavelet(5, :))/(TD_parameters.t(2) - TD_parameters.t(1)) ];
    synmodel.wavelet(5, :) = [ 0 diff(synmodel.wavelet(5, :))/(TD_parameters.t(2) - TD_parameters.t(1)) ];
    synmodel.wavelet(5, :) = [ 0 diff(synmodel.wavelet(5, :))/(TD_parameters.t(2) - TD_parameters.t(1)) ];
    synmodel.wavelet(5, :) = (synmodel.wavelet(5, :)/rms(synmodel.wavelet(5, :)));

    synmodel.polarization = linspace(0, 7*pi/8, nevt)';

    %per station parameters
    %synmodel.dt                = rand(nsta,2)*TD_parameters.max_dt;%ones(20, 1);% abs(randn(n,1)*TD_parameters.dt_std);
    %synmodel.fast_dir          = rand(nsta,2)*2*pi;%linspace(0, pi/2, 20)';%

    if TD_parameters.prior_information.use

        synmodel.dt(:, 1) = randn(nsta,1)*TD_parameters.prior_information.dt_std(1) + TD_parameters.prior_information.dt(1);
        synmodel.dt(:, 2) = randn(nsta,1)*TD_parameters.prior_information.dt_std(2) + TD_parameters.prior_information.dt(2);
        synmodel.fast_dir(:,1) = randn(nsta,1)*TD_parameters.prior_information.phi_std(1) + TD_parameters.prior_information.phi(1);%in radians
        synmodel.fast_dir(:,2) = randn(nsta,1)*TD_parameters.prior_information.phi_std(2) + TD_parameters.prior_information.phi(2);%in radians
        synmodel   = make_AB(synmodel);
        synmodel   = make_dtphi(synmodel);%enforce a pi/2 wrap
        synmodel.fast_dir_rotation(:,1) = randn(nsta,1)*TD_parameters.prior_information.rot_std(1) + TD_parameters.prior_information.rot(1);%in radians
        synmodel.fast_dir_rotation(:,2) = randn(nsta,1)*TD_parameters.prior_information.rot_std(2) + TD_parameters.prior_information.rot(2);%in radians

    else

        synmodel.A = randn(nsta,TD_parameters.max_layers);
        synmodel.B = randn(nsta,TD_parameters.max_layers);
        synmodel.fast_dir_rotation = randn(nsta,TD_parameters.max_layers)*TD_parameters.rotation_std;%in radians
        synmodel   = make_dtphi(synmodel);

    end

    synmodel.sta_or            = zeros(nsta,1);

    %per station-evt pairs
    synmodel.amp      = zeros(nevt, nsta);%log
    synmodel.shift    = randn(nevt, nsta)*TD_parameters.Delaystd;
    synmodel.dtS      = randn(nevt, nsta)*TD_parameters.dtS(2) + TD_parameters.dtS(1);%ones(20, 1);% abs(randn(n,1)*TD_parameters.dt_std);
    synmodel.sig      = rand(nevt, nsta)*0.5 + -3;

    %global parameters
    if TD_parameters.covariences 

        synmodel.r = TD_parameters.r_range(1) + randn()*TD_parameters.r_range(2);%rand(n,1)*TD_parameters.max_dt;%log(0.9);%
        synmodel.f = TD_parameters.f_range(1) + randn()*TD_parameters.f_range(2);%rand(n,1)*TD_parameters.max_dt;%log(0.5);

    else
        
        synmodel.r = TD_parameters.r_range(1);% + randn()*TD_parameters.r_range(2);%rand(n,1)*TD_parameters.max_dt;%log(0.9);%
        synmodel.f = TD_parameters.f_range(1);% + randn()*TD_parameters.f_range(2);%rand(n,1)*TD_parameters.max_dt;%log(0.5);

    end

    TD_parameters.polarization = synmodel.polarization;
    synmodel.polarization0     = synmodel.polarization;
    synmodel.polarization(:)   = 0;%no anomalies

    C = zeros(length(TD_parameters.t), length(TD_parameters.t));
    
    for i = 1:length(TD_parameters.t)

        for j = 1:length(TD_parameters.t)

            dt     = abs(TD_parameters.t(j) - TD_parameters.t(i));
            C(i,j) = exp(-(exp(synmodel.r))*dt)*cos(2*pi*exp(synmodel.f)*dt);

        end

    end
    
    window = tukeywin(length(TD_parameters.t), 0.01);
    
    prelogf = log([1; TD_parameters.f(2:end)]/(400*2*pi));

    for k = 1:nevt

        for kk = 1:nsta    

            [N, E] = apply_operator_multi(synmodel, TD_parameters.f, k, kk);
            
            N = new_waveform(N, synmodel.shift(k, kk), synmodel.amp(k, kk), synmodel.dtS(k, kk), TD_parameters.f, prelogf);
            E = new_waveform(E, synmodel.shift(k, kk), synmodel.amp(k, kk), synmodel.dtS(k, kk), TD_parameters.f, prelogf);

            noiseN = chol(C)*randn(size(TD_parameters.t))*exp(synmodel.sig(k, kk)).*window; 
            noiseE = chol(C)*randn(size(TD_parameters.t))*exp(synmodel.sig(k, kk)).*window; 

            n = N + noiseN;
            e = E + noiseE;

            norm = rms([ n; e]);
    
            snr = 1/(mean(abs(hilbert([ noiseN; noiseE]))) + 2*std(abs(hilbert([ noiseN; noiseE]))));%assumes ampllitude is 1

            allWfs(k,kk).north = n;%/norm;
            allWfs(k,kk).east  = e;%/norm;
            
            allWfs(k, kk).latitude      = 0;
            allWfs(k, kk).longitude     = 0;
            allWfs(k, kk).station       = [ 'Syn' num2str(kk) ];
            allWfs(k, kk).snr           = snr;
            allWfs(k, kk).orientation_k = 1./TD_parameters.orientation_std^2;

        end
        
    end

end