function model = center_wavelets(model, allWfs, TD_parameters)
    
    f = 2*pi*TD_parameters.f; %needs to be negative, but why?
    logf = log([1; f(2:end)]/(400*2*pi));

    [nevt, nsta] = size(allWfs);

    %center the wavelets to a delay of zero
    %set the amplitude
    %apply the opposite of the mean t* operator

    mdelay = mean(model.shift, 2);
    %mamp   = mean(model.amp, 2);
    mdtS   = mean(model.dtS, 2);

    for k = 1:nevt

        model.wavelet(k, :) = new_waveform(model.wavelet(k, :)', ...
            mdelay(k), 0, mdtS(k), TD_parameters.f, logf);

    end

    model = vectorize_model(model, TD_parameters);
    model = evaluate(model, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

end