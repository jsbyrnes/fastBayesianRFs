%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Driving script for SCEC RFs
%staind given on command line
addpath('Custom_functions_RF')
addpath('CircStat2012a')
addpath('irisFetch')
addpath('FetchData')
addpath('deconvolution_code/')
% addpath('./Reflectivity/toolbox/')
% addpath('./Reflectivity/deconvolution_code')
% addpath('./fastnonlinearACD/')
%addpath('./ParallelFastNonLinearACD-master/')
%addpath('./fminlbfgs_version2c')
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Define the data to be used 
name     = '3J_localevents100Hz';
dataName = '3J_nodes';
%%%%%%%%%%%%%%%%%%%%%

load('SCEC2019.mat')
evtind     = 48964; 
window     = 500;

rf_win = [ -0.25 2 ];

%%%%%%%%%%%%%%%%%%%%%
%configure the search
Parameters                 = make_parameters(name, dataName);
Parameters.parallel        = false;
Parameters.get_errors      = false;
Parameters.total_time      = 100;%in s
Parameters.solver_printout = false;%show progression of solver or not. Different for different solvers. 
% Parameters.low_pass        = 0.25;
% Parameters.sample_rate     = 3*0.25;%larger by at least 2

Parameters.reset_rounds      = 50;
Parameters.max_layers        = 1;
Parameters.rotation_std      = 1e-2*pi/180;

Parameters.low_pass    = 100;
Parameters.sample_rate = 300;%larger by at least 2

Parameters.high_pass   = 1;%3/(Parameters.total_time);%in seconds

Parameters.vertical_names    = { 'DPZ' };
Parameters.north_names       = { 'DPN' };
Parameters.east_names        = { 'DPE' };
Parameters.max_gaussians     = 20; 
Parameters.use_orientations  = false;%include orientations or not in model vector. 
Parameters.use_polarization  = false;
Parameters.use_covarience    = true;
Parameters.wavelet           = true;

Parameters.r_range             = [ log(0.5), 1 ];%mean, std, log
Parameters.f_range             = [ log(20), 2 ];%mean, std, log

Parameters.sig         = 0.25;
Parameters.niterations = 5e4;
Parameters.burnin      = 2.5e4;
Parameters.saveint     = 250;
Parameters.printon     = 1e3;
Parameters.update_sig  = 1e2;
Parameters.nchains     = 30;
Parameters.solver      = 'rjMCMC';
%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%
%Script is designed to use command line arguments to guide what data gets
%loaded
%Test if the test is synthetic by looking for the words "Syn" in the first
%three letters
%%%%%%%%
%rng('shuffle')
Parameters.t  = (0:1/Parameters.sample_rate:(Parameters.total_time))';

%frequency vector for FFTs
fs = 1/(Parameters.t(2) - Parameters.t(1));
dFreq = fs/length(Parameters.t); %frequency spacing
fNyq=fs/2;   %Nyquist frequency

if ~exist('station')

    station = [];

end

%the next two steps build the frequency vector (freqs at which spectrum
%was calculated) postivie and negative
f=(0:length(Parameters.t)-1)'*dFreq;
f(f>fNyq)=f(f>fNyq)-fNyq*2;
%Parameters.f = f;

if ~any(~(Parameters.dataName(1:3) == 'Syn'))

    %SynName and sig are a command line arguments
    
    if exist(['./' Parameters.dataName '/' Parameters.dataName '.mat' ])

        oldrun                  = load(['./' Parameters.dataName '/' Parameters.dataName '.mat' ], 'model', 'allWfs', 'synmodel', 'Parameters');
        allWfs                  = oldrun.allWfs;
        synmodel                = oldrun.synmodel;
        Parameters.polarization = oldrun.Parameters.polarization;%to avoid overwriting everything

    else

        %validate the forward problem here!!!
        [ allWfs, synmodel, Parameters ] = load_data_syn(Parameters);

    end

    %[ allWfs, t ] = load_data(Parameters, E(10), S);

else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    disp([ 'Getting the data for ' Parameters.dataName ]);

    if exist(['./' Parameters.name '/' Parameters.name 'Data.mat' ])

        oldrun                     = load(['./' Parameters.name '/' Parameters.name 'Data.mat' ], 'allWfs', 'Parameters');
        allWfs                     = oldrun.allWfs;

    else

        load([ './FetchData/' Parameters.dataName '.mat' ]);

        allWfs = load_SCEC_eq(Parameters, S, startTime(evtind), window, lat(evtind), lon(evtind));
    
        lat = [allWfs(:).latitude];
        lon = [allWfs(:).longitude];
    
        mkdir(['./' Parameters.name '/'])
        save(['./' Parameters.name '/' Parameters.name 'Data.mat' ], 'allWfs', 'Parameters');

    end

end

%pre-save kappas and orientations
[nevt, nsta] = size(allWfs);

for k = 1:nsta

    x   = [allWfs(:, k).orientation_k];

    if isempty(x(~isnan(x)))
        Parameters.sta_err(k,1) = 2*pi;
    else
        Parameters.sta_err(k,1) = unique(x(~isnan(x)));
        Parameters.sta_err(k,1) = sqrt(1/Parameters.sta_err(k,1));
    end

end

%count up how many traces you have (missing data are nans)
n = 0;
for k = 1:numel(allWfs)

    if ~any(isnan(allWfs(k).north))

        n = n + 1;

    end

end

Parameters.n = n;

t = Parameters.t;
% 
% figure(1)
% clf
% hold on
% 
% scale = range([allWfs(:).delta]);
% 
% for k = 39:39%length(allWfs)
% 
%     if ~isnan(allWfs(k).Z)
% 
%         plot(t - allWfs(k).t0, 0.1*scale*allWfs(k).Z/max(abs(allWfs(k).Z( (t - allWfs(k).t0)>rf_win(1) & (t - allWfs(k).t0)<rf_win(2) ))) + allWfs(k).delta, 'k')
%         plot(t - allWfs(k).t0, 0.1*scale*allWfs(k).R/max(abs(allWfs(k).Z( (t - allWfs(k).t0)>rf_win(1) & (t - allWfs(k).t0)<rf_win(2) ))) + allWfs(k).delta, 'r')
% 
%     end
% 
% end
% xlim([ -1 2 ])
%ylim([min([allWfs(:).delta])-0.1*scale, max([allWfs(:).delta])+0.1*scale])
%title('Vertical component')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%  opt                  = optimoptions('fminunc','Display','iter','Algorithm','quasi-newton',...
%     'MaxIterations', 1e9, 'OptimalityTolerance', eps, 'StepTolerance', 1e-4, ...
%     'MaxFunctionEvaluations', 1e9, 'FiniteDifferenceStepSize', 1e-6, 'FiniteDifferenceType', 'central');
%opt = optimset('display', 'iter', 'MaxFunEvals', 1e9, 'MaxIter', 1e9, 'TolFun', 0.1, 'TolX', 1e-4);

%r = randn()*Parameters.r_range(2) + Parameters.r_range(1);%rand(n,1)*Parameters.max_dt;%log(0.9);%
%f = randn()*Parameters.f_range(2) + Parameters.f_range(1);%rand(n,1)*Parameters.max_dt;%log(0.5);
%sig = randn() - 2;
%[rf, t, md, itnum ] = IDRFjsb('P', allWfs(k).Z, allWfs(k).R, 1/Parameters.sample_rate, -5, -5, 1, r, f, -1, 400);

disp([ 'On station #' num2str(sta_ind) ])

allWfs = allWfs(sta_ind);

P = allWfs.Z( (t - allWfs.t0)>rf_win(1) & (t - allWfs.t0)<rf_win(2) );
D = allWfs.R( (t - allWfs.t0)>rf_win(1) & (t - allWfs.t0)<rf_win(2) );

tw = tukeywin(length(P), 0.2);
P = P.*tw;
D = D.*tw;

[rfp, rf_mttmp] = multitaper2rf_3component(P, D, zeros(size(D)), 1/Parameters.sample_rate, ...
    1, length(D), 2.5, 3, 'P', [ Parameters.high_pass Parameters.low_pass ], [ 1 length(D)]);
%rf_timemt = (0:length(rf_mt)-1)/Parameters.sample_rate - 1;

[rf_td, ~] = IDRF('P', P, D, 1/Parameters.sample_rate, -1, -50, Parameters.low_pass, 1e-3, 1e-3, 400);

allWfs.Z = P;
allWfs.R = D;
%Parameters.t = Parameters.t( (t - allWfs(k).t0)>rf_win(1) & (t - allWfs(k).t0)<rf_win(2) );
%Parameters.t = Parameters.t - Parameters.t(1);

[rf_rjtmp, rf_stdtmp, ~ ] = fit_waveforms_rjMCMC(allWfs, Parameters);        

rf_rj  = [ zeros(1, Parameters.sample_rate) rf_rjtmp ];
rf_std = [ zeros(1, Parameters.sample_rate) rf_stdtmp ];
rf_mt  = [ rf_mttmp; zeros(Parameters.sample_rate,1) ];

%note - takes non-dimensional parameters
%bicfunc = @(x) IDRFjsb('P', P, D, 1/Parameters.sample_rate, -1, 10, 10, x(1), x(2), x(3), 400, Parameters);
%vec = fminunc(bicfunc, zeros(3,1), opt);
%vec = fminsearch(bicfunc, zeros(3,1), opt);
%best_sol = diff_evol(bicfunc, 20, -5*ones(3,1), 5*ones(3,1), 1e-4);
%[~, rf_btd, rf_timeb, ~, ~] = bicfunc(vec);
%bicfunc = @(x) IDRFjsb('P', allWfs(k).Z, allWfs(k).R, 1/Parameters.sample_rate, -5, 20, 1, x(1:2), x(3:4), x(5), 400, Parameters);
%vec = fminunc(bicfunc, randn(5,1), opt);
   
%     figure(2)
%     rf_time = (1:length(rf_rj))/Parameters.sample_rate - 1;
%     hold on
%     plot(rf_time, rf_mt(:, k))
%     plot(rf_time, rf_td(:, k))
%     plot(rf_time, rf_rj, 'k', 'LineWidth', 2)
%     plot(rf_time, rf_rj + rf_std/2, 'k--', 'LineWidth', 2)
%     plot(rf_time, rf_rj - rf_std/2, 'k--', 'LineWidth', 2)
%     xlim([-0.1 0.5])
%     xlabel('Time, s, relative to P arrival')
%     ylabel('Ps/P')
%     legend('Multitaper deconvolution', 'Iterative time domain deconvolution', 'rjMCMC deconvolution')

%model  = fit_waveforms(allWfs, Parameters);
%model = fit_waveforms_rjMCMC(allWfs, Parameters);

disp('-> Saving results')

%remove large fields that can be easily generated
% msave = [];
% for k = 1:length(model)
% 
%     m = rmfield(model(k), 'Cinv');
% 
%     msave = [ msave; m ];
% 
% end
% model = msave;

mkdir(['./' Parameters.name '/'])
if ~any(~(Parameters.name(1:3) == 'Syn'))

    save([ './' Parameters.name '/' Parameters.name num2str(rand(), 4) '.' Parameters.solver ...
        '.mat' ], 'allWfs', 'Parameters', 'model', 'synmodel')

else

    save([ './' Parameters.name '/' Parameters.name allWfs.station '.' Parameters.solver ...
        '.mat' ])

end
