function Parameters = make_parameters(name, dataName)

%%%%%%
%Editting this script is not recommended. 
%This script is best used as a "default" settings for the algorithm
%To customize for a run, edit in the execution script after calling this
%scipt

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Parameters.nth_root            = 3;%when making a guess at the wavelets
    Parameters.sig_range           = [ -3 2 ];%log. Traces are normalized so 0 is ~max, goes a bit high.
    Parameters.r_range             = [ log(0.05), 0.25 ];%mean, std, log
    Parameters.f_range             = [ log(0.05), 0.25 ];%mean, std, log
    Parameters.dtS                 = [ 0.25, 1 ];%log normal
    Parameters.tSAstd              = 0.25;%not log!
    Parameters.dtstd               = 0.5;%std
    Parameters.Delaystd            = 1.5;
    Parameters.max_layers          = 2;
    Parameters.polarization_std    = 5*pi/180;
    Parameters.orientation_std     = 5*pi/180;%default value
    Parameters.rotation_std        = 10*pi/180;%30 degree std for layer rotation
    Parameters.h                   = 1e-3;
    
    %%%%inverse problem parameters
    Parameters.solver            = 'bfgs';%current optiond are bfgs, acd, mcmc, paramonte
    Parameters.parallel          = false;
    Parameters.solver_printout   = false;%show progression of solver or not. Different for different solvers. 
    Parameters.use_orientations  = true;%include orientations or not in model vector. 
    Parameters.use_polarization  = true;
    Parameters.get_errors        = true;
    Parameters.use_tSA           = false;
    Parameters.use_covarience    = true;%include covarience parameters in model vector. 
    Parameters.wavelet           = true;%disable for cross-covolutional misfit

    %bfgs parameters
    Parameters.reset_size        = [ -1 -2 ];%non-dimensionalized, log10
    Parameters.reset_iter        = 5;%# of attempts without improvement
    Parameters.reset_rounds      = 25;
    
    %hmc parameters
    Parameters.nchains    = 10;
    Parameters.burnin     = 200;
    Parameters.batch_size = 100;
    Parameters.batch_thin = 10;
    
    Parameters.prior_information.use     = false;

    %data parameters
    Parameters.total_time  = 100;%in s
    Parameters.low_pass    = 1/12;
    Parameters.sample_rate = 0.3;%larger by at least 2
    Parameters.high_pass   = 1/33;%3/(Parameters.total_time);%in seconds
    Parameters.north_names = { 'BH1', 'HH1', 'BHN', 'HHN' };
    Parameters.east_names  = { 'BH2', 'HH2', 'BHE', 'HHE' };
    Parameters.name        = name;
    Parameters.dataName    = dataName;
    
    %misc, don't use. For specality tests
    Parameters.cluster           = false;%not recommended
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end


%%%%%%example of prior information
%     Parameters.prior_information.dt      = [ 1 1 ];
%     Parameters.prior_information.dt_std  = [ 0.33 0.33 ];
%     Parameters.prior_information.phi     = [ -80 78 ]*pi/180;
%     Parameters.prior_information.phi_std = [ 15 5 ]*pi/180;
%     Parameters.prior_information.rot     = [ 5 0  ]*pi/180;
%     Parameters.prior_information.rot_std = [ 15 5 ]*pi/180;
