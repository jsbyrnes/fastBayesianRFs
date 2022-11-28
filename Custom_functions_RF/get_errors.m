function model = get_errors(model, allWfs, TD_parameters)

    disp('-> Getting the Hessian for uncertainities')

    model = fill_hessian(model, allWfs, TD_parameters);

    model.C = inv(-model.H);%covarience matrix

end

%     c   = [];
%     df  = @(x) lpdf(x, model);
% 
%     disp('     Drawing samples with HMC')
%     for q = 1:hmc_chains
% 
%         smp = hmcSampler(df, randn(size(model.vector)), 'StepSize', 0.005);%very fast, so don't over think the numbers
%         tmp = drawSamples(smp, 'Burnin', 500, 'NumSamples', 100, 'ThinSize', 25);
%         c   = [ c; tmp];
% 
%     end
% 
%     model.correlations = corrcoef(c);
% 
%     %save errors
%     for k = 1:model.nparam
% 
%         [index, field, scale] = map_index(k, model, allWfs, TD_parameters);
%         field                 = [ field '_sigma' ];
%         model.(field)(index)  = std(c(:, k))*scale;
% 
%     end
% 
%     %reshape fields
%     model.amp_sigma     = reshape(model.amp_sigma, size(model.amp));
%     model.shift_sigma   = reshape(model.shift_sigma, size(model.shift));
%     model.dtS_sigma     = reshape(model.dtS_sigma, size(model.dtS));
%     model.wavelet_sigma = reshape(model.wavelet_sigma, size(model.wavelet));
%     model.sig_sigma     = reshape(model.sig_sigma, size(model.sig));
%     model.A_sigma       = reshape(model.A_sigma, size(model.A));
%     model.B_sigma       = reshape(model.B_sigma, size(model.B));
%     model.fast_dir_rotation_sigma     = reshape(model.fast_dir_rotation_sigma, size(model.fast_dir_rotation));
