function model = make_starting_models_NM(model, TD_parameters, allWfs)

    t = TD_parameters.t;

    [nevt, nsta] = size(allWfs);

    %per evt parameters
    for k = 1:nevt

        if TD_parameters.wavelet

            model.wavelet = model.wavelet + randn(size(model.wavelet))/10;

        else

            model.wavelet = TD_parameters.syn_wavelet;

        end
            
        model.polarization(k,1) = model.polarization(k,1) + randn()*TD_parameters.polarization_std/10;

    end

    %model.sta_or            = model.sta_or + randn(size(model.sta_or))/100;%very smooth in this case

    %per station-evt pairs
    model.sig      = model.sig + randn(size(model.sig))*TD_parameters.sig_range(2)/10;

    %global parameters
    model.r = model.r + TD_parameters.r_range(2)/10;
    model.f = model.f + TD_parameters.f_range(2)/10;

    model = build_C(model, t);
    model = evaluate(model, allWfs, TD_parameters, 1:nevt, 1:nsta, []);

    model = vectorize_model(model, TD_parameters);

end
