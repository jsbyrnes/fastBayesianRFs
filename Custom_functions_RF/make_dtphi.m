function model = make_dtphi(model)

    %vanderbeek and faccenda 2021
    model.dt       = sqrt(model.A.^2 + model.B.^2);
    model.fast_dir = atan(model.B./(model.dt + model.A));

end