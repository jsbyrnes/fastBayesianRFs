function model = apply_update(model, field, t)

    if strcmp(field, 'A') || strcmp(field, 'B')

        model = make_dtphi(model);

    elseif strcmp(field, 'dt') || strcmp(field, 'fast_dir')
    
        model = make_AB(model);
        model = make_dtphi(model);%does the wrapping correctly for fast_dir

    elseif strcmp(field, 'r') || strcmp(field, 'f')

        %model.r = min([ model.r log(0.5) ]);
        %model.f = min([ model.f log(1/(t(2) - t(1))) ]);%sampling rate the min

        model = build_C(model, t);

    end

end
