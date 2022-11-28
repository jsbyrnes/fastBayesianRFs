function model = make_AB(model)

    %vanderbeek and faccenda 2021

    n1 = sqrt(model.dt).*cos(model.fast_dir);
    n2 = sqrt(model.dt).*sin(model.fast_dir);

    model.A = (n1.^2 - n2.^2);
    model.B = 2*n1.*n2;

end