function [AIC, RF, RF_Time, itnum, D_cur] = IDRFjsb(Phase,P,D,dt,t_for,t_max,gauss_t,r,f,sig,itmax, Parameters)
% Iterative Deconvolution and Receiver-Function Estimation in time domain
% Version used in Junlin Hua's Sp RF + CCP stacking code, with some cleaning.

%Modified by JSB for bayesian analysis

% Phase is 'Sp' or 'Ps'
% P for parent phase (S in Sp case), D for daughter phase.
% dt: time between samples
% gauss_t: 1 std for the Gaussian convolved
% accept_mis1, accept_mis2 and itmax for stop criteria of the loop
% t_for,t_trun: remove amplitude with t > 0 in the Sp case 

    f   = exp(f*Parameters.f_range(2) + Parameters.f_range(1));%rand(n,1)*Parameters.max_dt;%log(0.5));
    r   = exp(r*Parameters.r_range(2) + Parameters.r_range(1));%rand(n,1)*Parameters.max_dt;%log(0.9));
    sig = exp(sig*Parameters.sig_range(2) + Parameters.sig_range(1));

    C = zeros(length(P), length(P));

    for i = 1:length(P)
    
        for j = 1:length(P)
    
            %dt     = abs(i - j)*dt;
            %C(i,j) = exp(-r*abs(i - j)*dt);%This form is highly stable when inverting
    
            for k = 1:length(f)

                C(i,j) = C(i,j) + exp(-r(k)*(abs(i - j))*dt)*cos(2*pi*f(k)*abs(i - j)*dt);
                %C(i,j) = C(i,j) + exp(-0.5*r(k)*(abs(i - j)*dt)^2)*cos(2*pi*f(k)*abs(i - j)*dt);
                %C(i,j) = C(i,j) + (cos(2*pi*f(k)*abs(i - j)*dt)/(1 + ((i-j)^2)/r^2 ));
                %C(i,j) = C(i,j) + (cos(2*pi*f(k)*abs(i - j)*dt).*sech(r(k)*(i-j)));

            end

        end
    
    end

    Cinv = inv(C);
    ld   = logdet(C, 'chol');

    md = @(v1, v2, t) (v1 - circshift(v2, t))'*(Cinv/sig^2)*(v1 - circshift(v2, t));
    %md = @(v1, v2, t) (v1 - circshift(v2, t))'*Cinv*(v1 - circshift(v2, t));

    RF_tmp  = zeros(length(P)*2-1,1);
    itnum   = 0;
    D_cur   = D;

    AIC     = 1e9;
    AIC_old = 1e10;

    index_vec = (-1*(floor(length(P)) - 1):(floor(length(P)) - 1));
    %[~,index_vec] = xcorr(D_cur,P);
    RF_Time = index_vec*dt;

    %while (misfit_old-misfit>accept_mis1*misfit || misfit>accept_mis2*misfit_ref) && itnum <= itmax
    while itnum <= itmax && AIC <= AIC_old 
        
        itnum = itnum+1;
        
%        if itnum == 1

%            ind = find(index_vec==0);

%        else

        %[amp_corr,~] = xcorr(D_cur,P);
        %auto_corr = xcorr(P);
        %[~,ind] = max(abs(amp_corr));
        %amp_rf = amp_corr(ind)/auto_corr((length(t_corr)+1)/2);
        mdiff = arrayfun( @ (x) md(D_cur, P, x), index_vec);
        %mdiffn = arrayfun( @ (x) md(D_cur, -1*P, x), -1*(floor(length(P)/2) - 1):(floor(length(P)/2) - 1));

        %test for both positive and negative - max misfit is negative amp
        mdiff    = (mdiff - mean(mdiff)).*tukeywin(length(mdiff), 0.5)';
        [~, ind] = max(abs(mdiff));
        %ind      = ind + (round(length(P)/2));

%        end

        %need to remake function each time to reload the new RF and ind
        added_RF = @(a) RF_tmp + return_delta(length(RF_tmp), ind, a);
        amp_fit = @(a) md(D, conv(P,added_RF(a),'same'),0);

        %misfit_old      = misfit;
        [ amp, misfit ] = fminbnd(amp_fit, -2, 2);
        RF_tmp          = added_RF(amp);

        %AIC is -2*llh + k*ln(n)
        %with 2 parameters per spike


        %RF_tmp(ind) = RF_tmp(ind)+amp_rf;
        D_sub = conv(P,RF_tmp,'same');
        D_cur = D - D_sub;

        llh = -0.5*log(2*pi)*length(D) - ld/2 - length(D)*log(sig) ...
               - misfit/2;

        AIC_old = AIC;
        AIC     = -2*llh + 2*(2*itnum + 3*length(f));

        %misfit_old = misfit;
        %misfit = sqrt(sum(D_cur.^2))/misfit_ref;

    end

    RF = RF_tmp;

    if strcmp(Phase,'Sp')
        RF(RF_Time>t_for)=0;
        RF = RF(RF_Time<=t_for);
        RF_Time = RF_Time(RF_Time<=t_for);
    else
        RF(RF_Time<t_for)=0;
        RF = RF(RF_Time>=t_for);
        RF_Time = RF_Time(RF_Time>=t_for);
    end

    RF      = RF(RF_Time<t_max);
    RF_Time = RF_Time(RF_Time<t_max);

    if gauss_t ~= 0
%         gauss_sig = gauss_t/dt;
%         x = linspace(-gauss_sig*4,gauss_sig*4,gauss_sig*8);
%         Gauss_win = exp(-x.^2/(4*gauss_sig^2));
%         RF = conv(RF,Gauss_win,'same');

        Nfft                   = length(RF);         % number  of points in fft = n of samples
        dF                     = 1/(dt*Nfft);       % frequency interval
        Nyq                    = (1/(2*dt));        % Nyquist frequency
        freqVals               = (0:(Nfft-1))*dF;         % this gives us frequencies with the correct spacing but going all the way to the sampling frequency
        freqVals(freqVals>Nyq) = freqVals(freqVals>Nyq)-(Nyq*2);
        %freqVals               = abs(freqVals);
        
        Faux = gauss_t;%number is half width in Hz
        A    = exp(-((freqVals.^2)/(2*(Faux^2))));
        A    = A'; %for dimensional consistency with inTr.data
        
        %Apply to the traces and revert to time domain
        amp  = max(RF);%max is a realiable spike
        RF   = real(ifft(fft(RF).*A));
        RF   = amp*RF/max(RF);

    end
    
end

function delta = return_delta(N, ind, a)

    delta = zeros(N,1);

    delta(ind) = a;

end