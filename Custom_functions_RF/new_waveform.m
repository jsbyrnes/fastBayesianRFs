function new = new_waveform(time_series, shift, A, tS, f, prelogf)

    tS = exp(tS);

    Amp = exp(-f*tS/2);                 % amplitude spectrum
    %%%%need to normalize amp
    Dt  = -(tS/pi)*prelogf; % dummy F=0, OK on next line
    Dt  = Dt - min(Dt);
    Dph = f.*(Dt);                   % delta phase shift spectrum
    % attn opperator in freq with a delay in time
    new = exp(A)*real(ifft(fft(time_series).*Amp.*exp(1i*Dph).*exp(-1i*2*pi*f*shift)));%operator for both the fast and slow directions

end
