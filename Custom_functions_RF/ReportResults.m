mkdir(['./' Parameters.name])
Parameters.t  = (0:1/Parameters.sample_rate:(Parameters.total_time))';
t = Parameters.t;

srcest = 0;

[~, indb] = max([model.lpst]);
m = model(indb);
m = build_C(m, Parameters.t);

[nevt, nsta] = size(allWfs);

if nsta > 1 && Parameters.wavelet

    %%%%%%%%%%%%%%%%%%%%%
    %plot traces by event
    for n = 1:nevt
    
        h = figure('visible', 'off');
        subplot(121)
        hold on
        subplot(122)
        hold on
    
        for k = 1:nsta
        
            if ~isnan(allWfs(n, k).north)
    
                [~, N, E, dN, dE] = apply_model(m, allWfs, Parameters, n, k);
                    
                amp = max([ N; E; dN; dE]);
    
                subplot(121)
                hold on
                plot(t, 0.5*dN/amp + k, 'k')
                plot(t, 0.5*N/amp + k, 'r')
        
                subplot(122)
                hold on
                plot(t, 0.5*dE/amp + k, 'k')
                plot(t, 0.5*E/amp + k, 'r')
    
            end
            
        end
    
        subplot(121)
        hold on
        title('North data (black) vs ML synthetic (red)')
        ylim([0 (nsta + 1)])
        subplot(122)
        hold on
        title('East data (black) vs ML synthetic (red)')
        ylim([0 (nsta + 1)])
    
        saveas(h, ['./' Parameters.name '/Event#' num2str(n) '.pdf'])
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%
    %plot traces by station
    
    for n = 1:nsta
    
        names = {allWfs(:, n).station};
        names(cellfun(@(names) any(isnan(names)),names)) = [];
    
        if isempty(names)
    
            continue %no data
    
        end
    
        staname = names{1};%will never be more than one
    
        h = figure('visible', 'off');
        subplot(121)
        hold on
        subplot(122)
        hold on
        
        for k = 1:nevt
        
            if ~isnan(allWfs(k,n).north)
                    
                [~, N, E, dN, dE] = apply_model(m, allWfs, Parameters, k, n);
                
                amp = max([ N; E; dN; dE]);
    
                subplot(121)
                hold on
                plot(t, 0.5*dN/amp + k, 'k')
                plot(t, 0.5*N/amp + k, 'r')
        
                subplot(122)
                hold on
                plot(t, 0.5*dE/amp + k, 'k')
                plot(t, 0.5*E/amp + k, 'r')
    
            end
            
        end
    
        subplot(121)
        hold on
        title('North data (black) vs ML synthetic (red)')
        ylim([0 (nevt + 1)])
        subplot(122)
        hold on
        title('East data (black) vs ML synthetic (red)')
        ylim([0 (nevt + 1)])
    
        saveas(h, ['./' Parameters.name '/' staname '_traces.pdf'])
    
    end
    close all

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlayers = length(model(1).dt(1,:));

%now make estimates of the parameters
for k = 1:numel(model)

    if k == 1

        for j = 1:nlayers

            fast_dir{j} = model(k).fast_dir(:, j);
            dt{j}       = model(k).dt(:, j);
            %A{j}        = model(k).A(:, j);
            %B{j}        = model(k).B(:, j);
            rotation{j} = model(k).fast_dir_rotation(:, j);
            tSA{j}      = model(k).tSA(:, j);

        end

        if nlayers > 1

            mphi    = circ_mean(2*model(k).fast_dir, [], 2)/2;
            diffphi = circ_dist(2*model(k).fast_dir(:, 2), 2*model(k).fast_dir(:, 1))/2;

        end

        r            = model(k).r;
        f            = model(k).f;

        if nsta > 1 && Parameters.wavelet

            polarization = model(k).polarization';
            sta_or       = model(k).sta_or';

        end

    else

        for j = 1:nlayers

            fast_dir{j} = [ fast_dir{j} model(k).fast_dir(:, j) ];
            dt{j}       = [ dt{j}       model(k).dt(:, j)       ];
            %A{j}        = [ A{j}        model(k).A(:, j)        ];
            %B{j}        = [ B{j}        model(k).B(:, j)        ];
            rotation{j} = [ rotation{j} model(k).fast_dir_rotation(:, j) ];
            tSA{j}      = [ tSA{j} model(k).tSA(:, j)           ];

        end

        if nlayers > 1

            mphi    = [ mphi circ_mean(2*model(k).fast_dir, [], 2)/2 ];
            diffphi = [ diffphi circ_dist(2*model(k).fast_dir(:, 2), 2*model(k).fast_dir(:, 1))/2 ];

        end

        r            = [ r model(k).r ];
        f            = [ f model(k).f ];

        if nsta > 1 && Parameters.wavelet

            polarization = [ polarization; model(k).polarization' ];
            sta_or       = [ sta_or; model(k).sta_or' ];

        end

    end

end

if Parameters.use_polarization && nsta > 1 && Parameters.wavelet

    h = figure('Visible', 'off');
    
    for k = 1:nevt
    
        mean_p(k) = circ_mean( (polarization(:, k)))*180/pi;
        std_p(k)  = circ_std(  (polarization(:, k)))*180/pi;
    
    end
    
    subplot(211)
    errorbar(mean_p, std_p, 'ko')
    xlabel('Event  #')
    ylabel('Polarization errors')
    subplot(212)
    plot(mean_p, std_p, 'ko')
    set(gca, 'XAxisLocation', 'origin')
    xlabel('Mean polarization error')
    ylabel('Mean polarization error standard deviation')
    saveas(h, ['./' Parameters.name '/PolarizationErrors.pdf'])

end

if Parameters.use_orientations && nsta > 1 && Parameters.wavelet
    
    h = figure('Visible', 'off');
    for k = 1:nsta
    
        mean_s(k) = circ_mean(sta_or(:, k))*180/pi;
        std_s(k)  = circ_std(sta_or(:, k))*180/pi;
    
    end
    
    subplot(211)
    errorbar(mean_s, std_s, 'ko')
    xlabel('Station #')
    ylabel('Orientation errors')
    subplot(212)
    plot(mean_s, std_s, 'ko')
    set(gca, 'XAxisLocation', 'origin')
    xlabel('Mean station orientation error')
    ylabel('Mean station orientation error standard deviation')
    saveas(h, ['./' Parameters.name '/OrientationErrors.pdf'])

end

if Parameters.use_tSA 

    for k = 1:nlayers

        tSA_l = tSA{k};

        tSA_l = mean(tSA_l);

    end

end

close all

if Parameters.get_errors

    %plot all the lpst results
    h = figure('Visible', 'off');
    subplot(121)
    lpst = [model(:).lpst];
    histogram(lpst(2:end))
    yl = ylim;
    hold on
    plot( [lpst(1) lpst(1)], [ 0 1e9 ], 'k-', 'LineWidth', 3)
    ylim(yl)
    xlim([min(lpst) (max(lpst) + 5)])
    xlabel('Log-posterior')
    title('HMC ensemble (histogram) and optimum (line)')
    axis square
    
    subplot(122)
    rv = log10(exp([model(:).r]));
    fv = log10(exp([model(:).f]));
    dr = linspace(min(rv), max(rv), 100);
    df = linspace(min(fv), max(fv), 100);
    [R,F] = meshgrid(dr, df);
    xrf   = [R(:), F(:)];
    ftp   = ksdensity([ rv; fv]', xrf, 'Function', 'pdf');%this is a little slow, but not too bad
    
    contourf(R, F, reshape(ftp, size(R)), 20), colormap(flipud(hot))
    xlabel('log_1_0(r)')
    ylabel('log_1_0(f)')
    title('Covarience parameters')
    axis square
    
    saveas(h, ['./' Parameters.name '/ModelQualities.pdf'])

    %use this to look at individual results
    %figure, k = 20; plot(dt(k, :), mod(fast_dir(k, :), pi), 'k.'), xlim([0, 4]), ylim([0, pi])
    
    %now test how many layers you have per station
    
    if Parameters.cluster
    
        nplot = 1;
    
    else
    
        nplot = nsta;
    
    end
    
    for k = 1:nplot
    
        h = figure('Visible', 'off');
    
        if Parameters.cluster
    
            staname = 'ClusteredData';
    
        else
    
            names = {allWfs(:, k).station};
            names(cellfun(@(names) any(isnan(names)),names)) = [];
        
            if isempty(names)
        
                rownames{k} = 'No Data';
                continue %no data
        
            end
        
            staname = names{1};%will never be more than one
            rownames{k} = staname;
    
        end
    
        ddt  = -2:0.025:6;%make this large to avoid problems with "reflection"
        dphi = (-pi/2):0.025:(pi/2);
        dr   = -40:2.5:40;
    
        [T1, P1]   = meshgrid(ddt, dphi);
        [T2, R2]   = meshgrid(ddt, dr);
        [P3, R3]   = meshgrid(dphi, dr);
        [PP1, PP2] = meshgrid(dphi, dphi);
        [TT1, TT2] = meshgrid(ddt, ddt);
        [RR1, RR2] = meshgrid(dr, dr);
    
        xi1 = [T1(:),  P1(:)];
        xi2 = [T2(:),  R2(:)];
        xi3 = [P3(:),  R3(:)];
        xi4 = [PP1(:), PP2(:)];
        xi5 = [TT1(:), TT2(:)];
        xi6 = [RR1(:), RR2(:)];
    
        ddtt  = -2:0.025:8;%make this large to avoid problems with "reflection"
        dmphi = (-pi/2):0.025:(pi/2);
        ddphi = (-pi):0.025:(pi);
        [tt1, mp]   = meshgrid(ddtt, dmphi);
        [tt2, dp]   = meshgrid(ddtt, ddphi);
    
        xi7 = [tt1(:), mp(:)];
        xi8 = [tt2(:), dp(:)];
    
        if nlayers == 1
            
            %get data for station, layer
            %Plots for 1 layer
            %%%%%%%%%%%%%%%%%%%%%%%
            %time against direction
            dtl       = dt{1};
            fast_dirl = fast_dir{1};
        
            dtl       = dtl(k, :);
            fast_dirl = fast_dirl(k, :);
    
            ftp = ksdensity([ dtl; fast_dirl]', xi1, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            ftp = reshape(ftp, size(T1));
        
            ftp = ftp/sum(ftp(:)*((pi/length(dphi))*(8/length(ddt))));%normalize the area to 1
        
            pdf_levels = linspace(0, max(ftp(:)), 1000);
        
            for n = 1:length(pdf_levels)
        
                ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
        
                %area under the curve
                area(n) = sum(ftp(ind)*(pi/length(dphi))*(8/length(ddt)));
        
            end
        
            levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
        
            %[~, ind] = max(ftp(:));
            %[iy, ix] = ind2sub(size(ftp), ind);
        
            subplot(1, 3, 1)
            contourf(T1, P1*180/pi, ftp, 7), colormap(flipud(hot))
            hold on
            contour(T1, P1*180/pi, ftp, [ levels(1) 1e9 ], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
            %contour(T1, P1*180/pi, ftp, [ levels(2) 1e9 ], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
            %plot(ddt(ix), dphi(iy)*180/pi, 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
            if exist('synmodel')
        
                plot(synmodel.dt(k, 1), (synmodel.fast_dir(k, 1))*180/pi, 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
            end
        
            xlim([0 4])
            xlabel('Splitting time, s')
            ylabel('Fast direction')
            axis square
            %%%%%%%%%%%%%%%%%%%%%%
            %time against rotation
            dtl       = dt{1};
            rotationl = rotation{1}*180/pi;
            
            dtl       = dtl(k, :);
            rotationl = rotationl(k, :);
    
            ftp = ksdensity([dtl; rotationl]', xi2, 'Function', 'pdf');%this is a little slow, but not too bad
            %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            ftp = reshape(ftp, size(R2));
        
            ftp = ftp/sum(ftp(:)*((200/length(dr))*(8/length(ddt))));%normalize the area to 1
        
            pdf_levels = linspace(0, max(ftp(:)), 1000);
        
            for n = 1:length(pdf_levels)
        
                ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
        
                %area under the curve
                area(n) = sum(ftp(ind)*(200/length(dr))*(8/length(ddt)));
        
            end
        
            levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
        
            %[~, ind] = max(ftp(:));
            %[iy, ix] = ind2sub(size(ftp), ind);
        
            subplot(1, 3, 2)
            contourf(T2, R2, ftp, 7), colormap(flipud(hot))
            hold on
            contour(T2, R2, ftp, [ levels(1) 1e9], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
            %contour(T2, R2, ftp, [ levels(2) 1e9], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
            %plot(ddt(ix), dr(iy), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
            if exist('synmodel')
        
                plot(synmodel.dt(k, 1), synmodel.fast_dir_rotation(k, 1)*180/pi, 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
            end
        
            xlabel('Splitting time, s')
            ylabel('Rotation, degrees')
            xlim([0 4])
            axis square
        
            %%%%%%%%%%%%%%%%%%%%%%
            %direction against rotation
            fast_dirl = fast_dir{1};
            rotationl = rotation{1}*180/pi;
        
            fast_dirl = fast_dirl(k, :);
            rotationl = rotationl(k, :);
    
            ftp = ksdensity([ fast_dirl; rotationl]', xi3, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            ftp = reshape(ftp, size(P3));
        
            ftp = ftp/sum(ftp(:)*((pi/length(dphi))*(200/length(dr))));%normalize the area to 1
        
            pdf_levels = linspace(0, max(ftp(:)), 1000);
        
            for n = 1:length(pdf_levels)
        
                ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
        
                %area under the curve
                area(n) = sum(ftp(ind)*(pi/length(dphi))*(200/length(dr)));
        
            end
        
            levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
        
            %[~, ind] = max(ftp(:));
            %[iy, ix] = ind2sub(size(ftp), ind);
        
            subplot(1, 3, 3)
            contourf(P3*180/pi, R3, ftp, 7), colormap(flipud(hot))
            hold on
            contour(P3*180/pi, R3, ftp, [ levels(1) 1e9], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
            %contour(P3*180/pi, R3, ftp, [ levels(2) 1e9], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
            %plot(dphi(ix), dr(iy), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
            if exist('synmodel')
        
                plot((synmodel.fast_dir(k, 1))*180/pi, synmodel.fast_dir_rotation(k, 1)*180/pi, 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
            end
        
            xlabel('Fast direction')
            ylabel('Rotation, degrees')
            axis square
        
            sgtitle([ 'For ' staname ])
            saveas(h, ['./' Parameters.name '/' staname '_1layer_PDF.pdf'])
            close all
    
            %%%%%%%%%%%%
            %layer 1
            dt1 = dt{1};
            fd  = fast_dir{1};
            rot = rotation{1};
    
            %optimum, mean/std of HMC, 95% internal for non-circular quantities
            layer1_dt_opt(k)  = model(1).dt(k, 1);
            layer1_dt_mean(k) = mean(dt1(k, :));
            layer1_dt_std(k)  = std(dt1(k, :));
            layer1_dt_25(k)   = quantile(dt1(k, :),0.025);
            layer1_dt_975(k)  = quantile(dt1(k, :),0.975);
    
            layer1_fd_opt(k)  = model(1).fast_dir(k, 1)*180/pi;
            layer1_fd_mean(k) = (180/pi)*circ_mean(2*fd(k, :)')/2;%note correction for 2 thetaness
            layer1_fd_std(k)  = (180/pi)*circ_std(2*fd(k, :)')/2;
            
            layer1_rot_opt(k)  = model(1).fast_dir_rotation(1)*180/pi;
            layer1_rot_mean(k) = mean(rot(k, :))*180/pi;
            layer1_rot_std(k)  = std(rot(k, :))*180/pi;
            layer1_rot_25(k)   = quantile(rot(k, :),0.025)*180/pi;
            layer1_rot_975(k)  = quantile(rot(k, :),0.975)*180/pi;
    
        elseif nlayers == 2 
    
            h = figure('Visible', 'off');
            for j = 1:nlayers
                
                %plots for two layer model
                %%%%%%%%%%%%%%%%%%%%%%%
                %time against direction
                dtl       = dt{j};
                fast_dirl = fast_dir{j};
            
                dtl       = dtl(k, :);
                fast_dirl = fast_dirl(k, :);
    
                ftp = ksdensity([ dtl; fast_dirl]', xi1, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
                %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
                ftp = reshape(ftp, size(T1));
                
                ftp = ftp/sum(ftp(:)*((pi/length(dphi))*(8/length(ddt))));%normalize the area to 1
                
                pdf_levels = linspace(0, max(ftp(:)), 1000);
                
                for n = 1:length(pdf_levels)
                
                    ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
                
                    %area under the curve
                    area(n) = sum(ftp(ind)*(pi/length(dphi))*(8/length(ddt)));
                
                end
                
                levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
                
                %[~, ind] = max(ftp(:));
                %[iy, ix] = ind2sub(size(ftp), ind);
                    
                subplot(nlayers+1, 3, (j-1)*3 + 1 )
                contourf(T1, P1*180/pi, ftp, 7), colormap(flipud(hot))
                hold on
                contour(T1, P1*180/pi, ftp, [ levels(1) 1e9], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
                %contour(T1, P1*180/pi, ftp, [ levels(2) 1e9], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
                %plot(ddt(ix), dphi(iy)*180/pi, 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
                if exist('synmodel')
        
                    plot(synmodel.dt(k, j), (synmodel.fast_dir(k, j))*180/pi, 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
                end
        
                xlim([0 4])
                title([ 'Layer #' num2str(j) ])
                xlabel('Splitting time, s')
                ylabel('Fast direction')
        
                %%%%%%%%%%%%%%%%%%%%%%
                %time against rotation
                dtl       = dt{j};
                rotationl = rotation{j}*180/pi;
    
                dtl       = dtl(k, :);
                rotationl = rotationl(k, :);
    
                ftp = ksdensity([ dtl; rotationl]', xi2, 'Function', 'pdf');%this is a little slow, but not too bad
                %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
                ftp = reshape(ftp, size(R2));
                
                ftp = ftp/sum(ftp(:)*((200/length(dr))*(8/length(ddt))));%normalize the area to 1
                
                pdf_levels = linspace(0, max(ftp(:)), 1000);
                
                for n = 1:length(pdf_levels)
                
                    ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
                
                    %area under the curve
                    area(n) = sum(ftp(ind)*(200/length(dr))*(8/length(ddt)));
                
                end
                
                levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
                
                %[~, ind] = max(ftp(:));
                %[iy, ix] = ind2sub(size(ftp), ind);
                    
                subplot(nlayers+1, 3, (j-1)*3 + 2 )
                contourf(T2, R2, ftp, 7), colormap(flipud(hot))
                hold on
                contour(T2, R2, ftp, [ levels(1) 1e9], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
                %contour(T2, R2, ftp, [ levels(2) 1e9], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
                %plot(ddt(ix), dr(iy), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
                if exist('synmodel')
        
                    plot(synmodel.dt(k, j), synmodel.fast_dir_rotation(k, j)*180/pi, 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
                end
        
                title([ 'Layer #' num2str(j) ])
                xlabel('Splitting time, s')
                ylabel('Rotation, degrees')
                xlim([0 4])
        
                %%%%%%%%%%%%%%%%%%%%%%
                %direction against rotation
                fast_dirl =  fast_dir{j};
                rotationl = rotation{j}*180/pi;
        
                fast_dirl = fast_dirl(k, :);
                rotationl = rotationl(k, :);
                
                ftp = ksdensity([ fast_dirl; rotationl]', xi3, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
                %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
                ftp = reshape(ftp, size(P3));
                
                ftp = ftp/sum(ftp(:)*((pi/length(dphi))*(200/length(dr))));%normalize the area to 1
                
                pdf_levels = linspace(0, max(ftp(:)), 1000);
                
                for n = 1:length(pdf_levels)
                
                    ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
                
                    %area under the curve
                    area(n) = sum(ftp(ind)*(pi/length(dphi))*(200/length(dr)));
                
                end
                
                levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
                
                %[~, ind] = max(ftp(:));
                %[iy, ix] = ind2sub(size(ftp), ind);
                    
                subplot(nlayers+1, 3, (j-1)*3 + 3 )
                contourf(P3*180/pi, R3, ftp, 7), colormap(flipud(hot))
                hold on
                contour(P3*180/pi, R3, ftp, [ levels(1) 1e9], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
                %contour(P3*180/pi, R3, ftp, [ levels(2) 1e9], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
                %plot(dphi(ix)*180/pi, dr(iy), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
                if exist('synmodel')
        
                    plot((synmodel.fast_dir(k, j))*180/pi, synmodel.fast_dir_rotation(k, j)*180/pi, 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
        
                end
        
                title([ 'Layer #' num2str(j) ])
                xlabel('Fast direction')
                ylabel('Rotation, degrees')
        
            end
        
            %%%%%%%%%%%%%%%%%%%%%%
            %time against time
            dt1 =  dt{1};
            dt2 =  dt{2};
        
            dt1 = dt1(k, :);
            dt2 = dt2(k, :);
            
            ftp = ksdensity([ dt1; dt2]', xi5, 'Function', 'pdf');%this is a little slow, but not too bad
            %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            ftp = reshape(ftp, size(TT1));
        
            ftp = ftp/sum(ftp(:)*((8/length(ddt))^2));%normalize the area to 1
        
            pdf_levels = linspace(0, max(ftp(:)), 1000);
        
            for n = 1:length(pdf_levels)
        
                ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
        
                %area under the curve
                area(n) = sum(ftp(ind)*(((8/length(ddt))^2)));
        
            end
        
            levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
        
            %[~, ind] = max(ftp(:));
            %[iy, ix] = ind2sub(size(ftp), ind);
        
            subplot(nlayers+1, 3, 7 )
            contourf(TT1, TT2, ftp, 7), colormap(flipud(hot))
            hold on
            contour(TT1, TT2, ftp, [ levels(1) 1e9], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
            %contour(TT1, TT2, ftp, [ levels(2) 1e9], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
            %plot(ddt(ix), ddt(iy), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
    
            if exist('synmodel')
    
                plot((synmodel.dt(k, 1)), synmodel.dt(k, 2), 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
    
            end
    
            xlabel('Splitting Time, layer 1')
            ylabel('Splitting Time, layer 2')
            xlim([0 4]); ylim([0 4]);
        
            %%%%%%%%%%%%%%%%%%%%%%
            %direction against direction
            fast_dir1 =  fast_dir{1};
            fast_dir2 =  fast_dir{2};
        
            fast_dir1 = fast_dir1(k, :);
            fast_dir2 = fast_dir2(k, :);    
    
            ftp = ksdensity([ fast_dir1; fast_dir2]', xi4, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            ftp = reshape(ftp, size(PP1));
        
            ftp = ftp/sum(ftp(:)*((pi/length(dphi))^2));%normalize the area to 1
        
            pdf_levels = linspace(0, max(ftp(:)), 1000);
        
            for n = 1:length(pdf_levels)
        
                ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
        
                %area under the curve
                area(n) = sum(ftp(ind)*((pi/length(dphi))^2));
        
            end
        
            levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
        
            %[~, ind] = max(ftp(:));
            %[iy, ix] = ind2sub(size(ftp), ind);
        
            subplot(nlayers+1, 3, 8 )
            contourf(PP1*180/pi, PP2*180/pi, ftp, 7), colormap(flipud(hot))
            hold on
            contour(PP1*180/pi, PP2*180/pi, ftp, [ levels(1) 1e9], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
            %contour(PP1*180/pi, PP2*180/pi, ftp, [ levels(2) 1e9], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
            %plot(dphi(ix)*180/pi, dphi(iy)*180/pi, 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
    
            if exist('synmodel')
    
                plot((synmodel.fast_dir(k, 1))*180/pi, (synmodel.fast_dir(k, 2))*180/pi, 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
    
            end
    
            title([ 'Layer #' num2str(j) ])
            xlabel('Fast dir, layer 1')
            ylabel('Fast dir, layer 2')
            
            %%%%%%%%%%%%%%%%%%%%%%
            %rotation against rotation
            dr1 =  rotation{1}*180/pi;
            dr2 =  rotation{2}*180/pi;
        
            dr1 = dr1(k, :);
            dr2 = dr2(k, :);    
    
            ftp = ksdensity([ dr1; dr2]', xi6, 'Function', 'pdf');%this is a little slow, but not too bad
            %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            ftp = reshape(ftp, size(RR1));
        
            ftp = ftp/sum(ftp(:)*((200/length(dr))^2));%normalize the area to 1
        
            pdf_levels = linspace(0, max(ftp(:)), 1000);
        
            for n = 1:length(pdf_levels)
        
                ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
        
                %area under the curve
                area(n) = sum(ftp(ind)*(((200/length(dr))^2)));
        
            end
        
            levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
        
            %[~, ind] = max(ftp(:));
            %[iy, ix] = ind2sub(size(ftp), ind);
        
            subplot(nlayers+1, 3, 9)
            contourf(RR1, RR2, ftp, 7), colormap(flipud(hot))
            hold on
            contour(RR1, RR2, ftp, [ levels(1) 1e9], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
            %contour(RR1, RR2, ftp, [ levels(2) 1e9], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
            %plot(dr(ix), dr(iy), 'x', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
    
            if exist('synmodel')
    
                plot((synmodel.fast_dir_rotation(k, 1))*180/pi, (synmodel.fast_dir_rotation(k, 2))*180/pi, 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
    
            end
    
            title([ 'Layer #' num2str(j) ])
            xlabel('Rotation, layer 1')
            ylabel('Rotation, layer 2')
        
            sgtitle([ '2 layer PDFs for ' staname ])
            saveas(h, ['./' Parameters.name '/' staname '_2layers_PDF.pdf'])
            close all
    
            %%%%%%%%%%%%%%%
            %now make seperate plots for the "collected" parameters
    
            h = figure('Visible', 'off');
    
            %%%%%%%%%%%%%%%%%%%%%%
            %toal time against mean directions
            subplot(211)
            dtt   =  dt{1} + dt{2};
        
            dtt   = dtt(k, :);
            mphik = mphi(k, :);    
    
            ftp = ksdensity([ dtt; mphik]', xi7, 'Function', 'pdf');%this is a little slow, but not too bad
            %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            ftp = reshape(ftp, size(tt1));
        
            ftp = ftp/sum(ftp(:)*(pi/length(dmphi))*(10/length(ddtt)));%normalize the area to 1
        
            pdf_levels = linspace(0, max(ftp(:)), 1000);
        
            for n = 1:length(pdf_levels)
        
                ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
        
                %area under the curve
                area(n) = sum(ftp(ind)*(pi/length(dmphi))*(10/length(ddtt)));
        
            end
        
            levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
            
            contourf(tt1, mp*180/pi, ftp, 7), colormap(flipud(hot))
            hold on
            contour(tt1, mp*180/pi, ftp, [ levels(1) 1e9], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
            %contour(tt1, mp*180/pi, ftp, [ levels(2) 1e9], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
    
            if exist('synmodel')
    
                plot(synmodel.dt(k, 1) + synmodel.dt(k, 2), 0.5*circ_mean(synmodel.fast_dir(k, 1)*2, synmodel.fast_dir(k, 2)*2)*180/pi, 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
    
            end
    
            xlabel('Total Splitting Time')
            ylabel('Mean fast direction')
            xlim([0 8]); ylim([-90 90]); grid on
    
            %toal time against change in directions
            subplot(212)
            dtt   =  dt{1} + dt{2};
        
            dtt   = dtt(k, :);
            dphik = diffphi(k, :);    
    
            ftp = ksdensity([ dtt; dphik]', xi8, 'Function', 'pdf');%this is a little slow, but not too bad
            %f = ksdensity([ dt(k, :); wrapToPi(fast_dir(k, :)*2)/2]', xi, 'Function', 'pdf', 'BoundaryCorrection','reflection');%this is a little slow, but not too bad
            ftp = reshape(ftp, size(tt2));
        
            ftp = ftp/sum(ftp(:)*(2*pi/length(ddphi))*(10/length(ddtt)));%normalize the area to 1
        
            pdf_levels = linspace(0, max(ftp(:)), 1000);
        
            for n = 1:length(pdf_levels)
        
                ind = ftp(:)>pdf_levels(n);%how much of the pdf is above this level?
        
                %area under the curve
                area(n) = sum(ftp(ind)*(2*pi/length(ddphi))*(10/length(ddtt)));
        
            end
        
            levels = interp1(area + cumsum(1e-8*ones(size(area))), pdf_levels, [ 0.68, 0.95]);
            
            contourf(tt2, dp*180/pi, ftp, 7), colormap(flipud(hot))
            hold on
            contour(tt2, dp*180/pi, ftp, [ levels(1) 1e9], 'LineWidth', 3, 'LineColor', 'k', 'LineStyle', '--')
            %contour(tt2, dp*180/pi, ftp, [ levels(2) 1e9], 'LineWidth', 1, 'LineColor', 'k', 'LineStyle', '--')
    
            if exist('synmodel')
    
                plot(synmodel.dt(k, 1) + synmodel.dt(k, 2), 0.5*circ_dist(2*synmodel.fast_dir(k, 2), 2*synmodel.fast_dir(k, 1))*180/pi, 'o', 'MarkerSize', 20, 'LineWidth', 5, 'Color', [0.5 0.5 0.5])
    
            end
    
            xlabel('Total Splitting Time')
            ylabel('Differences in fast directions')
            xlim([0 8]); ylim([-180 180]); grid on
    
            sgtitle([ 'Collected PDFs for ' staname ])
            saveas(h, ['./' Parameters.name '/' staname '_Collected_PDF.pdf'])
            close all
    
            %%%%%%%%%%%
            %save the results in a table
    
            %layer 1
            dta  = dt{1};
            fda  = fast_dir{1};
            rota = rotation{1};
    
            %optimum, mean/std of HMC, 95% internal for non-circular quantities
            layer1_dt_opt(k)  = model(1).dt(k, 1);
            layer1_dt_mean(k) = mean(dta(k, :));
            layer1_dt_std(k)  = std(dta(k, :));
            layer1_dt_25(k)   = quantile(dta(k, :),0.025);
            layer1_dt_975(k)  = quantile(dta(k, :),0.975);
    
            layer1_fd_opt(k)  = model(1).fast_dir(k, 1)*180/pi;
            layer1_fd_mean(k) = (180/pi)*circ_mean(2*fda(k, :)')/2;%note correction for 2 thetaness
            layer1_fd_std(k)  = (180/pi)*circ_std(2*fda(k, :)')/2;
            
            layer1_rot_opt(k)  = model(1).fast_dir_rotation(k, 1)*180/pi;
            layer1_rot_mean(k) = mean(rota(k, :))*180/pi;
            layer1_rot_std(k)  = std(rota(k, :))*180/pi;
            layer1_rot_25(k)   = quantile(rota(k, :),0.025)*180/pi;
            layer1_rot_975(k)  = quantile(rota(k, :),0.975)*180/pi;
    
            %layer 2
            dta  = dt{2};
            fda  = fast_dir{2};
            rota = rotation{2};
    
            %optimum, mean/std of HMC, 95% internal for non-circular quantities
            layer2_dt_opt(k)  = model(1).dt(k, 2);
            layer2_dt_mean(k) = mean(dta(k, :));
            layer2_dt_std(k)  = std(dta(k, :));
            layer2_dt_25(k)   = quantile(dta(k, :),0.025);
            layer2_dt_975(k)  = quantile(dta(k, :),0.975);
    
            layer2_fd_opt(k)  = model(1).fast_dir(k, 2)*180/pi;
            layer2_fd_mean(k) = (180/pi)*circ_mean(2*fda(k, :)')/2;%note correction for 2 thetaness
            layer2_fd_std(k)  = (180/pi)*circ_std(2*fda(k, :)')/2;
            
            layer2_rot_opt(k)  = model(1).fast_dir_rotation(k, 2)*180/pi;
            layer2_rot_mean(k) = mean(rota(k, :))*180/pi;
            layer2_rot_std(k)  = std(rota(k, :))*180/pi;
            layer2_rot_25(k)   = quantile(rota(k, :),0.025)*180/pi;
            layer2_rot_975(k)  = quantile(rota(k, :),0.975)*180/pi;
    
        end
    
    end
    
    T = table(layer1_dt_opt', layer1_dt_mean', layer1_dt_std', layer1_dt_25', layer1_dt_975', ...
        layer1_fd_opt', layer1_fd_mean', layer1_fd_std', ...
        layer1_rot_opt', layer1_rot_mean', layer1_rot_std', layer1_rot_25', layer1_rot_975', ...
        'VariableNames', { 'Splitting Time, s - optimum', 'Splitting Time, s - HMC mean', 'Splitting Time, s - HMC std', ...
        'Splitting Time, s - 95% lower bound', 'Splitting Time, s - 95% upper bound', 'Fast direction, ° - optimum', ...
        'Fast direction, ° - HMC mean', 'Fast direction, ° - HMC std', 'Rotation, ° - optimum', 'Rotation, ° - HMC mean', ...
        'Rotation, ° - HMC std','Rotation, ° - 95% lower bound','Rotation, ° - 95% upper bound' }, ...
        'RowNames', rownames');
    
    writetable(T, ['./' Parameters.name '/Layer1' ], 'FileType', 'spreadsheet', 'WriteRowNames', true)
    
    if nlayers == 2
    
        T = table(layer2_dt_opt', layer2_dt_mean', layer2_dt_std', layer2_dt_25', layer2_dt_975', ...
            layer2_fd_opt', layer2_fd_mean', layer2_fd_std', ...
            layer2_rot_opt', layer2_rot_mean', layer2_rot_std', layer2_rot_25', layer2_rot_975', ...
            'VariableNames', { 'Splitting Time, s - optimum', 'Splitting Time, s - HMC mean', 'Splitting Time, s - HMC std', ...
            'Splitting Time, s - 95% lower bound', 'Splitting Time, s - 95% upper bound', 'Fast direction, ° - optimum', ...
            'Fast direction, ° - HMC mean', 'Fast direction, ° - HMC std', 'Rotation, ° - optimum', 'Rotation, ° - HMC mean', ...
            'Rotation, ° - HMC std','Rotation, ° - 95% lower bound','Rotation, ° - 95% upper bound' }, ...
            'RowNames', rownames');
    
        writetable(T, ['./' Parameters.name '/Layer2' ], 'FileType', 'spreadsheet', 'WriteRowNames', true)
    
    end

end
