function [ allWfs, synmodel, TD_parameters ] = make_NoMelt_syn(TD_parameters)

    %arbitary plotting, shift them around
    %shift1      = 27;%in s. Works for high f
    shift1      = 0;%in s. Works for high f
    
    GorW = 1;%1 is geographic (north - east), 2 is wave (radial transverse)
    
    baz_vector  = 0:5:180;
    
    dbaz = baz_vector(2) - baz_vector(1);
        
    dz   = 2.5;
    vpvs = 1.76;
    
    fid = fopen('./Reflectivity/Data7621_bw1_05025.40.0.-103.0.csv');
    fgetl(fid);fgetl(fid);fgetl(fid);
    bids = fgetl(fid);
    fclose(fid);
    
    bids = strsplit(bids, ',');
    
    indM = str2num(bids{1}) + 2;
    indN = str2num(bids{2}) + 2;
    
    m = readmatrix('./Reflectivity/Data7621_bw1_05025.40.0.-103.0.csv', 'NumHeaderLines', 4);
    
    m(:, 2) = cumsum(m(:, 2));
    
    Ndepth = 0.5*(m(indN, 2) + m(indN-1, 2));
    
    %%%%%%%%%%%
    %make a new model in the format of the model structure
    depth = dz:dz:350;
    
    model.z    = depth;
    model.vp   = vpvs*interp1(m(:, 2), m(:, 1), depth, 'linear', 'extrap');
    
    % model.z    = [ 10, 75, 125, 175, 225, 230];
    % model.vp   = 8*ones(size(model.z));
    % model.vpvs = vpvs*ones(size(model.vp));
    
    %now fill out with zeros, not doing anisotropy here, then density
    model.theta = 90*ones(size(model.vp));
    model.phi   = 90*ones(size(model.vp));
    model.A     = zeros(size(model.vp));
    model.B     = zeros(size(model.vp));
    model.C     = zeros(size(model.vp));
    model.rho   = nafedrake_rho(model.vp); %g/cm^3
    
    % model.C   = [0, 0.1, 0.05, 0.15, 0.025, 0];
    % model.phi = [0,  10, 45,    -45, 20,    0 ];
    % 
    load('./Reflectivity/NoMelt_ANISO_med');
    
    G   = sqrt( Gs.^2 + GcGPa.^2 );
    psi = 0.5*atan(Gs./GcGPa)*180/pi;
        
    depthiso = depthiso(9:end) - depthiso(7);
    Vsv      = Vsv(9:end)/1000;
    
    Vsv(:) = 4.3; 
    
    %Load the NoMelt model
    model.vp   = interp1(depthiso + cumsum(0.001*ones(size(depthiso))), Vsv*vpvs, model.z, 'linear', 'extrap');
    model.rho  = nafedrake_rho(model.vp); %g/cm^3
    model.vpvs = vpvs*ones(size(model.vp));
    
    %G(depth>100) = 0;
    %%%%%%%%%%%%%%%%%
    %apply the actual NoMelt model
%     model.phi          = interp1(depth, psi, model.z, 'linear', 0) + 90;
%     model.C            = interp1(depth, G*1e9, model.z, 'linear', 0);%G in Pa
%     
%     %convert to peak to peak velocity
%     L = (1000*model.rho.*(1000*model.vp./model.vpvs).^2);
%     
%     model.C = 2*(sqrt((L + interp1(depth, G*1e9, model.z, 'linear', 0))./(1000*model.rho)) ...
%         - sqrt((L)./(1000*model.rho)))./sqrt((L)./(1000*model.rho));
%     
%     model.C(length(model.z)-3:end) = 0;%isotropic halfspace and then some

    %%%%%%%%%%%%%%%%
    %simpler two layer model
    model.C(:) = 0.00;
    % % model.C(model.z > 50 & model.z < 150) = linspace(0.1, 0.0, ...
    % %    length(model.phi(model.z > 50 & model.z < 150)));
    % %model.C(model.z > 50 & model.z < 100) = 0.1;
    % % model.C(model.z >= 100 & model.z < 150) = 0.05;
    % 
    model.C(model.z > 0 & model.z <= 50)   = 0.05;
    model.C(model.z > 50 & model.z <= 100) = 0.075;
    
    % model.phi(model.z > 0 & model.z <= 100) = wrapTo360(linspace(30, 210, ...
    %      length(model.phi(model.z > 0 & model.z <= 100))));
    
    model.phi(model.z > 0 & model.z <= 50) = wrapTo360(linspace(60, 45, ...
        length(model.phi(model.z > 0 & model.z <= 50))));
    
    model.phi(model.z > 50 & model.z <= 100) = wrapTo360(linspace(70, 100, ...
        length(model.phi(model.z > 50 & model.z <= 100))));

    %model.C(:) = 0;%isotropic halfspace and then some
    %layer thickness
    thk = diff([ 0 model.z]);
        
    for k = 1:length(baz_vector)
        
        %now set other parameters needed to make a receiver function
        baz    = baz_vector(k);
        phase  = 'SV';
        slow   = 0.01;
        dt     = .05;
        
        [P_comp, SV_comp, SH_comp] = anirec(phase, dt, slow, baz, model, 0);
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %make a source that looks like a teleseismic S wave in frequency
        Nfft                   = length(P_comp);         % number  of points in fft = n of samples
        dF                     = 1/(dt*Nfft);       % frequency interval
        Nyq                    = (1/(2*dt));        % Nyquist frequency
        freqVals               = (0:(Nfft-1))*dF;         % this gives us frequencies with the correct spacing but going all the way to the sampling frequency
        freqVals(freqVals>Nyq) = freqVals(freqVals>Nyq)-(Nyq*2);
        %freqVals               = abs(freqVals);
        
        t = (0:(length(P_comp) - 1))'*dt;

         w = (exp(-((t - mean(t) - 6).^2)/(2*3.^2)) ...
             + exp(-((t - mean(t) + 6).^2)/(2*3.^2)));
%         w = (exp(-((t - mean(t) - 2).^2)/(2*(0.25.^2))));
         w = [ 0; diff(w)/(t(2) - t(1)) ];
        A = w/rms(w);%source

%         Faux = source(2)/sqrt(log(2));%number is half width in Hz
%         A    = exp(-((freqVals-source(1))/Faux).^2);
%         A    = A'; %for dimensional consistency with inTr.data
        
%         %Apply to the traces and revert to time domain
%         P_comp  = real(ifft(fft(P_comp).*A));
%         SV_comp = real(ifft(fft(SV_comp).*A));
%         SH_comp = real(ifft(fft(-1*SH_comp).*A));
        
    %     SV_comp = lowpass(SV_comp, corner, 1/0.05);
    %     SH_comp = lowpass(SH_comp, corner, 1/0.05);
    
        %[~, ind] = max(sqrt(SV_comp.^2 + SH_comp.^2));
    
        SV_comp = real(ifft(fft(SV_comp).*fft(A)));
        SH_comp = real(ifft(fft(SH_comp).*fft(A)));

        SV_comp = circshift(SV_comp, -round(shift1/dt));
        SH_comp = circshift(-SH_comp, -round(shift1/dt));    
    
        SV_comp(t>TD_parameters.total_time) = [];
        SH_comp(t>TD_parameters.total_time) = [];
        t(t>TD_parameters.total_time) = [];

        Nr = SV_comp*cosd(baz) + SH_comp*sind(baz);
        Er = SV_comp*sind(baz) - SH_comp*cosd(baz);

        N(:, k) = interp1(linspace(0, TD_parameters.total_time, length(Nr)),...
            Nr, TD_parameters.t, 'pchip');
        E(:, k) = interp1(linspace(0, TD_parameters.total_time, length(Er)),...
            Er, TD_parameters.t, 'pchip');

        amp = rms([ N(:, k); E(:, k)]);

        N(:, k) = N(:, k)/amp;
        E(:, k) = E(:, k)/amp;

    end
    
    t = (0:(length(P_comp) - 1))'*dt;
    %A = circshift(A, -950);%for the high f wavelet
    A(t>TD_parameters.total_time) = [];
    TD_parameters.syn_wavelet = interp1(linspace(0, TD_parameters.total_time, length(A)),...
        A, TD_parameters.t, 'pchip');

    %global parameters
    %synmodel.r = -1;%randn()*TD_parameters.r_range(2) + TD_parameters.r_range(1);%rand(n,1)*TD_parameters.max_dt;%log(0.9);%
    %synmodel.f = log(1);%randn()*TD_parameters.f_range(2) + TD_parameters.f_range(1);%rand(n,1)*TD_parameters.max_dt;%log(0.5);
    synmodel.r = randn()*TD_parameters.r_range(2) + TD_parameters.r_range(1);%rand(n,1)*TD_parameters.max_dt;%log(0.9);%
    synmodel.f = randn()*TD_parameters.f_range(2) + TD_parameters.f_range(1);%rand(n,1)*TD_parameters.max_dt;%log(0.5);

    TD_parameters.polarization = baz_vector'*pi/180;

    C = zeros(length(TD_parameters.t), length(TD_parameters.t));
    
    for i = 1:length(TD_parameters.t)
        for j = 1:length(TD_parameters.t)
            dt     = abs(TD_parameters.t(j) - TD_parameters.t(i));
            C(i,j) = exp(-( exp(synmodel.r)/exp(synmodel.f))*dt)*cos(2*pi*exp(synmodel.f)*dt);
        end
    end
    
    window = tukeywin(length(TD_parameters.t), 0.01);

    for k = 1:length(baz_vector)

        noiseN = chol(C)*randn(size(TD_parameters.t))*exp(-2).*window;
        noiseE = chol(C)*randn(size(TD_parameters.t))*exp(-2).*window;

        n = N(:, k) + noiseN;
        e = E(:, k) + noiseE;

        norm = rms([ n; e]);

        snr = 1/(mean(abs(hilbert([ noiseN; noiseE]))) + 2*std(abs(hilbert([ noiseN; noiseE]))));%assumes ampllitude is 1

        allWfs(k,1).north = n;%/norm;
        allWfs(k,1).east  = e;%/norm;

        allWfs(k, 1).latitude      = 0;
        allWfs(k, 1).longitude     = 0;
        allWfs(k, 1).station       = [ 'Syn' num2str(1) ];
        allWfs(k, 1).snr           = snr;
        allWfs(k, 1).orientation_k = 1./TD_parameters.orientation_std^2;
        
    end

end