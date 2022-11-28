function model = make_starting_models(Parameters, allWfs)

    [nevt, nsta] = size(allWfs);

    model.t = ((1:length(allWfs(1).Z))/Parameters.sample_rate)';
    %per evt parameters

    model.n = 1;%randi([ 1 Parameters.max_gaussians ], [nsta, 1]);

    model.A{1} = 0.4;%randn(model.nz(k),1);
    model.t0   = 0.0;%randn(model.nz(k),1);
    model.w    = 2*log(model.t(2));%randn(model.nz(k),1);
    %model.R.dp   = 0;

    %model.T.A{1} = 0.4;%randn(model.nz(k),1);
    %model.T.t    = 0.0;%randn(model.nz(k),1);
    %model.T.w    = log(dt);%randn(model.nz(k),1);
    %model.T.dp   = 0;

    %per station-evt pairs
    model.sig      = ones(nevt, nsta);%*Parameters.sig_range(2) + Parameters.sig_range(1);%zeros(nevt, nsta)*diff(Parameters.sig_range) + Parameters.sig_range(1);
    %model.sigT      = ones(nevt, nsta);%*Parameters.sig_range(2) + Parameters.sig_range(1);%zeros(nevt, nsta)*diff(Parameters.sig_range) + Parameters.sig_range(1);

    %global parameters
    %start with just one
    model.r = randn()*Parameters.r_range(2) + Parameters.r_range(1);%rand(n,1)*Parameters.max_dt;%log(0.9);%
    model.f = randn()*Parameters.f_range(2) + Parameters.f_range(1);%rand(n,1)*Parameters.max_dt;%log(0.5);
    model.a = 1;

    model.temperature = exprnd(5) + 1;

    model.update_rf = true;

    model   = build_C(model);
    model   = evaluate(model, allWfs, Parameters);

end

% n = floor(length(Parameters.t)/3);
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

%     [nevt, nsta] = size(allWfs);
% 
%     model.t = Parameters.t;
%     dt      = model.t(2) - model.t(1);
%     minf    = 1/(model.t(end)*0.5);
%     maxf    = 0.5/dt;
%     %per evt parameters
% 
%     options                  = optimoptions('fminunc','Display', 'none', 'Algorithm','quasi-newton',...
%         'MaxIterations', 1e9, 'OptimalityTolerance', 1e-6, ...
%         'MaxFunctionEvaluations', 1e9);
% 
%     for k = 1:nevt
%     
%        if nsta > 1
% 
%            wavelet = zeros(length(Parameters.t),1);
%     
%             %stack the stack to the nth root
%             for j = 1:nsta
%         
%                 if ~any(isnan(allWfs(k,j).Z))
%     
%                     %wavelet = wavelet + sign(allWfs(k,j).Z).*abs(allWfs(k,j).Z).^(1/Parameters.nth_root);
%                     wavelet = wavelet + allWfs(k,j).Z;
%     
%                 end
%     
%             end
%     
%             %model.sources{k} = sign(wavelet).*abs(wavelet).^(Parameters.nth_root);
%             wavelet = wavelet/rms(wavelet);% + 0.25*randn(size(wavelet)));
%         
%             %shift back by 5 seconds
%             wavelet = circshift(wavelet, round(-5/dt));
%     
%             %now get an approximation to it (not super fast but not super slow)
%             a           = randn(2,1);
%             relativefit = sum((wavelet - harmonic_sum(a, length(model.t))).^2)/sum(wavelet.^2);
%     
%             while relativefit > 0.01
%     
%                 a = [ a; randn() ];
%     
%                 f = @(a) sum((wavelet - harmonic_sum(a, length(model.t))).^2);
%                 a = fminunc(f, a, options);
%     
%                 relativefit = f(a)/sum(wavelet.^2);
%     
%             end
%     
%             model.sources{k} = a;
% 
%        end
% 
%     end
% 
%     model.nrf = ones([nsta 1]);%randi([ 1 Parameters.max_gaussians ], [nsta, 1]);
% 
%     for k = 1:nsta
% 
%         model.az{k} = 1;%randn(model.nz(k),1);
%         model.Ar{k} = 0.4;%randn(model.nz(k),1);
%         model.Br{k} = 0.0;%randn(model.nz(k),1);
%         model.At{k} = 0.0;%randn(model.nz(k),1);
%         model.Bt{k} = 0.0;%randn(model.nz(k),1);
% 
%         model.f0{k} = minf;
% 
%         model.t0{k}  = 0;%allWfs(1,k).phase_times(1);%rand(model.nz(k),1)*model.t(end);
%         model.w{k}   = log(dt);%randn(model.nz(k),1);
% 
%         model.dp{k}  = 0;
% 
%         %model.at{k}  = 0;%randn(model.nz(k),1);
% 
%     end
% 
%     %per station-evt pairs
%     model.sigR      = ones(nevt, nsta);%*Parameters.sig_range(2) + Parameters.sig_range(1);%zeros(nevt, nsta)*diff(Parameters.sig_range) + Parameters.sig_range(1);
%     model.sigT      = ones(nevt, nsta);%*Parameters.sig_range(2) + Parameters.sig_range(1);%zeros(nevt, nsta)*diff(Parameters.sig_range) + Parameters.sig_range(1);
% 
%     %global parameters
%     model.r = randn()*Parameters.r_range(2) + Parameters.r_range(1);%rand(n,1)*Parameters.max_dt;%log(0.9);%
%     model.f = randn()*Parameters.f_range(2) + Parameters.f_range(1);%rand(n,1)*Parameters.max_dt;%log(0.5);
% 
%     model.temperature = 1 + exprnd(10);
% 
%     model.vector_source = true;
% 
%     model   = build_C(model, Parameters.t);
%     model   = evaluate(model, allWfs, Parameters, 1:nevt, 1:nsta);

