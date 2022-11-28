function [t0,t1,y0,y1,dt0,E,h_E] = envelopePick(time,s,tau,twb,twa,ramp,npts,niter,tol,varargin)
% ENVELOPEPICK--Autopicker for seismic arrival times. Arrival times are
% chosen by fitting a linear profile to the envelop of the waveform. This
% code is a modified version of the routine described in Shintaku et al.
% [2014].
%
% -B. VanderBeek
%
% Code trys to find the best-fit linear profile to the envelope of the
% seismic trace. The linear profile is defined by two horizontal lines and
% a connecting ramp:
%
%               (t1,y1).__________________
%                     /
% __________________./(t0,y0)
%
% Code iteratively inverts for the best-fit t0, y0, t1, y1 parameters where
% t0 is the arrival time.
%
% INPUT
% <Required>
%     time: Time vector corresponding to the seismic trace (s)
%        s: Seismic trace vector
%      tau: Initial pick time, i.e. t0 in above diagram. Setting tau equal
%           to range/velocity usually provides an adequet first guess (s).
%      twb: time window before ramp (s)
%      twa: time window after ramp (s)
%     npts: Width of the moving average filter applied to seismic trace
%           envelope in points
%     ramp: The initial ramp length, i.e. t1-t0 in the above diagram (s)
%    niter: The number of iterations (convergence usually happens in < 10
%           iterations)
% <Optional>
%   tf_verbose: If true, displays messages pertaining to solution. Default
%               is false.
%      tf_plot: If true, plots each solution. Default is false.
%
% OUTPUT
%    t0: Arrival time (s)
%    t1: Ending time of ramp (s)
%    y0: base of ramp (--)
%    y1: top of ramp (--)
%   dt0: Estimated standard error in arrival time (s)
%     E: Smoothed natural log of the waveform envelope
%
% NOTES
% [1] Suggested input values...
%      tau = range/velocity
%     npts = 10 pts
%      twa = 10 s 
%      twb = 10 s
%     ramp = 2 s
%    niter = 10
%
% REFERENCES
% [1] Shintaku, N., Forsyth, D. W., Hajewski, C. J., & Weeraratne, D. S. 
%     (2014). Pn anisotropy in Mesozoic western Pacific lithosphere. 
%     Journal of Geophysical Research: Solid Earth, 119(4), 3050-3063.
%

% Variable inputs
if length(varargin) == 1
    tf_verbose = varargin{1};
    tf_plot = false;
elseif length(varargin) == 2
    tf_verbose = varargin{1};
    tf_plot = varargin{2};
else
    tf_verbose = false;
    tf_plot = false;
end

% Define figure handles
if tf_plot
    lcolor = jet(niter);
    lcolor = [0.5*ones(1,3); lcolor];
    h_E    = figure; hold on;
    xlabel('time (s)');
    ylabel('envelope amplitude');
    title('Auto-pick Results');
else
    h_E = [];
end

% Force npts to be odd integer
npts = round(npts) + ~logical(rem(npts,2));

% Calculate and smooth the natural log of the envelope
E = abs(hilbert(s));
E = E - min(E) + 1; % log of zero is -Inf
E = log(E);
E = smooth(E,npts);
time = time(:); % force column vectors
E    = E(:);

%%%% Start auto-pick %%%%
% Trying to find the best-fit linear profile to the enveliope of the
% seismic trace

% Initial guess
%     (t1,y1).________
%           /
% ________./(t0,y0)
%

t0 = tau; y0 = mean(E) - std(E);
t1 = tau + ramp; y1 = mean(E) + std(E);
for kk = 1:niter
    try
        % Set pick window
        lla = (time >= t0-twb) & (time < t0);
        llb = (time >= t0) & (time < t1);
        llc = (time >= t1) & (time < t1+twa);
        % Define predicted data in each domain
        ya = y0*ones(size(time(lla)));
        yb = (((y1-y0)/(t1-t0))*(time(llb)-t0)) + y0;
        yc = y1*ones(size(time(llc)));
        % The predicted data
        yi = [ya(:);yb(:);yc(:)];
        
        % Illustrate inversion procedure
        if tf_plot
            figure(h_E);
            plot(time(logical(lla+llb+llc)),E(logical(lla+llb+llc)),'-k','linewidth',2); hold on;
            plot(time(logical(lla+llb+llc)),yi,'-','Color',lcolor(kk,:),'linewidth',2);
            xlim([t0-twb,t1+twa]);
            box on;
        end
        
        % Misfit
        dy = E(logical(lla+llb+llc))-yi;
        
        % Define partial derivatives in each domain
        % Pre-arrival
        dtdt0a = zeros(size(time(lla)));
        dtdy0a = ones(size(time(lla)));
        dtdt1a = zeros(size(time(lla)));
        dtdy1a = zeros(size(time(lla)));
        % Ramp up
        dtdt0b = (y1-y0)*(time(llb)-t1)./((t1-t0)^2);
        dtdy0b = ((t0-time(llb))/(t1-t0)) + 1;
        dtdt1b = (time(llb)-t0)*((y0-y1)/((t1-t0)^2));
        dtdy1b = (time(llb)-t0)/(t1-t0);
        % Post-arrival
        dtdt0c = zeros(size(time(llc)));
        dtdy0c = zeros(size(time(llc)));
        dtdt1c = zeros(size(time(llc)));
        dtdy1c = ones(size(time(llc)));
        
        % Make sure each domain is defined by at least 5 pts
        if (sum(lla) <= 5) || (sum(llb) <= 5) || (sum(llc) <= 5)
            error('Deficient');
        end
        
        % Define coefficient matrix and invert
        A = [[dtdt0a;dtdt0b;dtdt0c],[dtdy0a;dtdy0b;dtdy0c],[dtdt1a;dtdt1b;dtdt1c],[dtdy1a;dtdy1b;dtdy1c]];
        
        % Invert (a few options)
        % dx = A\dy; % direct
        % dx = lsqr(A,dy); % fast approximation of above
        [dx,stdx] = lscov(A,dy); % This gives estimated standard errors, though their magnitudes are ambiguous
        [~, msgid] = lastwarn; % Catch rank deficient warning
        
        % Update model parameters
        t0 = t0 + dx(1);
        y0 = y0 + dx(2);
        t1 = t1 + dx(3);
        y1 = y1 + dx(4);
        % Estimated error in pick time
        dt0 = stdx(1);
        
        % Forcing
        % This ensures the ramp does not reverse slope
        if t0 > t1
            t0 = t1-0.2;
            if tf_verbose
                display('Pick adjusted');
            end
        end
        
    catch
        if tf_verbose
            warning('No solution found!');
        end
        t0 = NaN;
        dt0 = NaN;
        msgid = '';
        break
    end
    
end
% Final solution plot
if ~isnan(t0)
    % Set pick window
    lla = (time >= t0-twb) & (time < t0);
    llb = (time >= t0) & (time < t1);
    llc = (time >= t1) & (time < t1+twa);
    % Define predicted data in each domain
    ya = y0*ones(size(time(lla)));
    yb = (((y1-y0)/(t1-t0))*(time(llb)-t0)) + y0;
    yc = y1*ones(size(time(llc)));
    % The predicted data
    yi = [ya(:);yb(:);yc(:)];
    % Illustrates how inversion proceeds
    if tf_plot
        figure(h_E);
        plot(time(logical(lla+llb+llc)),E(logical(lla+llb+llc)),'-k'); hold on;
        plot(time(logical(lla+llb+llc)),yi,'-','Color',lcolor(kk+1,:),'linewidth',2);
        xlim([t0-twb,t1+twa]);
        box on;
    end
end

% Check solution
if strcmp(msgid,'MATLAB:lscov:RankDefDesignMat')
    lastwarn('');
    t0 = NaN;
    dt0 = NaN;
end
if y1 < (y0 + tol*std(E(lla)));
    if tf_verbose
        warning(['Pick discarded because tolerance was not met. Ramp amplitude is within ',num2str(tol),'-standard deviations of noise level before arrival.']);
    end
    t0 = NaN;
    dt0 = NaN;
end
