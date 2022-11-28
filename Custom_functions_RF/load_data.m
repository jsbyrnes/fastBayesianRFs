function allWfs = load_data(Parameters, E, S, target_phase, pre)

    %some networks are right handed, but typical convections are left
    %handed. So for some listed here you have to flip stuff to make it work

    left_handed_networks = { 'YL' };

    stations_to_change = { 'G02B', 'G04B'};
    statiion_new_name  = { 'G02',  'G04' };

    if nargin == 4

        pre = 0.5*Parameters.total_time*ones(size(E));%centered

    end

    post = Parameters.total_time - pre;

    %first, how much padding to do to account for the 50% taper. 

    %check if there are orientations
    if exist([ './FetchData/' Parameters.dataName '_orientations.mat'])==2

        load([ './FetchData/' Parameters.dataName '_orientations.mat'], 'orientation', 'orientation_error', 'Station')

        rotate = true;

    else

        rotate = false;

    end

    channel_string = [];

    for k = 1:length(Parameters.vertical_names)
        channel_string = [ channel_string Parameters.vertical_names{k}, ','];
    end
    for k = 1:length(Parameters.north_names)
        channel_string = [ channel_string Parameters.north_names{k}, ','];
    end
    for k = 1:length(Parameters.east_names)
        channel_string = [ channel_string Parameters.east_names{k}, ','];
    end

    for evt = 1:length(E)
        
        t = Parameters.t;    
        
        for ks=1:length(S)
    
            %first, scale the array
            allWfs(evt, ks).north         = nan;
            allWfs(evt, ks).east          = nan;
            allWfs(evt, ks).Z             = nan;
            allWfs(evt, ks).R             = nan;
            allWfs(evt, ks).T             = nan;
            allWfs(evt, ks).latitude      = nan;
            allWfs(evt, ks).longitude     = nan;
            allWfs(evt, ks).snr           = nan;
            allWfs(evt, ks).station       = nan;
            allWfs(evt, ks).phase_list    = nan;
            allWfs(evt, ks).phase_times   = nan;
            allWfs(evt, ks).delta         = nan;
            allWfs(evt, ks).backazi       = nan;
            allWfs(evt, ks).orientation   = nan;
            allWfs(evt, ks).orientation_k = nan;
                        
            %first check to see if this station was active when the earthquake
            %happened
            
            sDate=datenum(S(ks).StartDate);
            if isempty(S(ks).EndDate); S(ks).EndDate='2500-01-01 00:00:00.000'; end
            eDate=datenum(S(ks).EndDate);
            eqDate=datenum(E(evt).PreferredTime);
            
            if ~(eqDate>sDate && eqDate<eDate) %if eq not between eDate and sDate
                %disp('No data for chosen event')
                continue
            end
            
            %calculate Delta
            [D,Az] = distance(S(ks).Latitude,S(ks).Longitude,E(evt).PreferredLatitude,E(evt).PreferredLongitude);
            Az     = Az - 180;
            phase_list = 'P,PP';%fixed, needed for identifying the actual exact phases
            urlstr     = sprintf('http://service.iris.edu/irisws/traveltime/1/query?distdeg=%f&evdepth=%f&noheader=true&mintimeonly=true&phases=%s',D,E(evt).PreferredDepth,phase_list);
                
            try
            
                tStr=urlread(urlstr);
        
                tCell=textscan(tStr,'%f %f %s %f %f %f %f %f %s %s');
            
            catch

                continue

            end

            phases= tCell{3};
            times = tCell{4};
            rayp  = tCell{5};

            if length(target_phase) == 1

                phaseTime = times(strcmpi(phases,target_phase{1}));
                rayp      = rayp(strcmpi(phases,target_phase{1}));

            else

                phaseTime = times(strcmpi(phases,target_phase{evt}));                
                rayp      = rayp(strcmpi(phases,target_phase{evt}));                

            end
    
            if isempty(phaseTime)

                continue %out of range
                
            end

            originTimeStr  = E(evt).PreferredTime;
            originTimeNum  = datenum(originTimeStr);
            
            startTime = originTimeNum + (phaseTime - pre(evt))/(24*60*60);            
            endTime   = originTimeNum + (phaseTime + post(evt))/(24*60*60);

            %add in padding. 
            pad_seconds = (pre(evt) + post(evt))*8;

            endTime   = endTime   + pad_seconds/(24*60*60);%accounts for a 50% taper in bandpassfilt
            startTime = startTime - pad_seconds/(24*60*60);           

            startTime = datestr(startTime,'yyyy-mm-dd HH:MM:SS.FFF');   
            endTime   = datestr(endTime,'yyyy-mm-dd HH:MM:SS.FFF');
        
            %to save later
            times = ((times/(24*60*60) + originTimeNum) - datenum(startTime))*(24*60*60) - pad_seconds;%sec relative to start of trace

            %fprintf('Trying event #%.0f station %s\n',ke,S(ks).StationCode)
            try
                
                %if strcmp(email, '')
                
                myTrace = irisFetch.Traces(S(ks).NetworkCode,...
                    S(ks).StationCode,'*',channel_string,startTime, endTime);
                        
    %             else
    %                
    %                 myTrace = irisFetch.Traces(S(ks).NetworkCode,...
    %                         S(ks).StationCode,'*',channel_string,startTime, endTime, ...
    %                         { email, password });
    %                 
    %             end
                    
            catch ME
    
                continue
                
            end
           
            %check the channels and keep highest sample rate with
            %dupilicates

            if length(myTrace)>3

                names = {myTrace.channel};

                %find the north channel (might be prerotation)
                vind = find(contains(names, Parameters.vertical_names));
                
                %find the north channel (might be prerotation)
                nind = find(contains(names, Parameters.north_names));

                %find the east channel (might be prerotation)
                eind = find(contains(names, Parameters.east_names));

                if length(vind) > 1

                    [~, ind] = max([myTrace(nind).sampleRate]);
                    vind = vind(ind);

                    [~, ind] = max([myTrace(nind).sampleRate]);
                    nind = nind(ind);

                    [~, ind] = max([myTrace(eind).sampleRate]);
                    eind = eind(ind);

                end

                myTrace = myTrace([nind eind]);

            end

            % Get rid of empty structures or missing data
            data  = {myTrace.data};
            tf_empty = cellfun('isempty',data);
            myTrace  = myTrace(~tf_empty);
    
            if length(myTrace)~=3
    
                continue
    
            end

            %check if an anomlously wrong amount of data was downloaded
            if abs(myTrace(1).startTime - datenum(startTime))*24*60*60 > 5

                continue

            end

            if abs(myTrace(1).endTime - datenum(endTime))*24*60*60 > 5

                continue

            end

            if abs(myTrace(2).startTime - datenum(startTime))*24*60*60 > 5

                continue

            end

            if abs(myTrace(2).endTime - datenum(endTime))*24*60*60 > 5

                continue

            end

            for q = 1:length(myTrace)
            
                %demean/detrend
                myTrace(q).data = detrend(myTrace(q).data);

            end

%             %not great, but there's a lot of padding for a reason
%             if length(myTrace(1).data) ~= length(myTrace(2).data)
% 
%                 if length(myTrace(1).data) > length(myTrace(2).data)
% 
%                     myTrace(2).data = [ myTrace(2).data; zeros(length(myTrace(1).data) - length(myTrace(2).data), 1) ];
% 
%                 elseif length(myTrace(1).data) < length(myTrace(2).data)
% 
%                     myTrace(1).data = [ myTrace(1).data; zeros(length(myTrace(2).data) - length(myTrace(1).data), 1) ];
% 
%                 end
% 
%             end

            try

                %causal response correction
                %myTrace = wfRemInstResp(myTrace);%filters

            catch

                error('Instrument response failed')

            end

            if ~isempty(myTrace)
                                                       
                %convert to second, divide by 4, and convert to index
                padding_ind = round(myTrace(1).sampleRate*pad_seconds);%convert to seconds and divide by 4

                if min([ myTrace(:).sampleCount ]) <= padding_ind

                    continue

                end

                %upper limits    
                myTrace = wfButterworth(myTrace, [ Parameters.high_pass Parameters.low_pass ]);

                for q = 1:length(myTrace)

                    %remove the padding
                    myTrace(q).data(1:padding_ind)                             = [];
                    myTrace(q).data((length(myTrace(q).data)-padding_ind):end) = [];

                end
    
                myTrace = wfResample_jsb(myTrace, Parameters);

                for q = 1:length(myTrace)
    
                    if length(t) > length(myTrace(q).data)
    
                        %not my favorite solution but almost always just
                        %one sample
                        myTrace(q).data = [ myTrace(q).data; myTrace(q).data(end)*ones(length(t) - length(myTrace(q).data), 1) ];
    
                    elseif length(t) < length(myTrace(q).data)
    
                        myTrace(q).data = myTrace(q).data(1:length(t));
    
                    end
    
                    %check if the channels are consistent at all
                    rms_chan(q) = rms(myTrace(q).data);

                end
            
%                 if rms_chan(1)/rms_chan(2) > 10 || rms_chan(2)/rms_chan(1) > 10
% 
%                     continue %bad data on a channel, remove. Relative amplitudes matter a lot.
% 
%                 end
                
                names = {myTrace.channel};
                
                data = [];

                %find the north channel (might be prerotation)
                vind = find(contains(names, Parameters.vertical_names),1);

                %find the north channel (might be prerotation)
                nind = find(contains(names, Parameters.north_names),1);

                %find the east channel (might be prerotation)
                eind = find(contains(names, Parameters.east_names),1);
                
                if ~isempty(vind) && ~isempty(nind) && ~isempty(eind)
    
                    %find the one that matches
                    %put in the right index                   

                    allWfs(evt, ks).Z     = myTrace(vind).data;

                    if contains(myTrace(1).network, left_handed_networks)

                        allWfs(evt, ks).north     = myTrace(eind).data;
                        allWfs(evt, ks).east      = myTrace(nind).data;

                    else

                        allWfs(evt, ks).north     = myTrace(nind).data;
                        allWfs(evt, ks).east      = myTrace(eind).data;

                    end

                    allWfs(evt, ks).latitude  = myTrace(eind).latitude;
                    allWfs(evt, ks).longitude = myTrace(eind).longitude;
                    allWfs(evt, ks).station   = myTrace(eind).station;

                    allWfs(evt, ks).phase_list  = phases;
                    allWfs(evt, ks).phase_times = times;%in secs into trace
                    allWfs(evt, ks).delta       = D;%in secs into trace
                    allWfs(evt, ks).backazi     = Az;%in secs into trace
                    allWfs(evt, ks).rayp        = rayp;%s/deg
                    
                    data = [myTrace(vind).data myTrace(nind).data myTrace(eind).data];%for normalization

                    if wrapTo360(myTrace(eind).azimuth - myTrace(nind).azimuth) == 90

                        dN = allWfs(evt,ks).north;
                        dE = allWfs(evt,ks).east;

                        allWfs(evt,ks).north = cosd(myTrace(nind).azimuth)*dN - sind(myTrace(nind).azimuth)*dE;
                        allWfs(evt,ks).east  = sind(myTrace(nind).azimuth)*dN + cosd(myTrace(nind).azimuth)*dE;

                    end

                else

                    continue

                end
    
                if rotate && ((myTrace(eind).azimuth - myTrace(nind).azimuth) ~= 90)

                    dN    = allWfs(evt,ks).north;
                    dE    = allWfs(evt,ks).east;
                    theta = orientation(strcmp(Station, allWfs(evt, ks).station));

                    if (isempty(theta) || isnan(theta))

                        allWfs(evt,ks).orientation = nan;
                        allWfs(evt,ks).orientation_k = 1/(100)^2;%effectively a uniform when put into von mises. Data will be hard to use.

                    else

                        %not actually north and east yet
                        allWfs(evt,ks).north = cosd(theta)*dN - sind(theta)*dE;
                        allWfs(evt,ks).east  = sind(theta)*dN + cosd(theta)*dE;

                        allWfs(evt,ks).orientation = theta*pi/180;

                    end

                    theta_std = orientation_error(strcmp(Station, allWfs(evt, ks).station))*pi/180;

                    if isempty(theta_std) || isnan(theta_std)

                        allWfs(evt,ks).orientation_k = 1/(100)^2;

                    else

                        theta_std                    = max([ Parameters.orientation_std, theta_std]);%dont let it be zero
                        allWfs(evt,ks).orientation_k = 1/(theta_std)^2;%note that its not a "standard deviation" anymore, but inverse varience

                    end

                else

                    if wrapTo360(myTrace(eind).azimuth - myTrace(nind).azimuth) == 90

                        allWfs(evt,ks).orientation   = 0;%orientated station
                        allWfs(evt,ks).orientation_k = 1/(Parameters.orientation_std)^2;%...maybe

                    else %ocean bottom seismometers are both listed as zero

                        allWfs(evt,ks).orientation   = 0;%number doesn't matter
                        allWfs(evt,ks).orientation_k = 1/(100)^2;%basically a uniform distribution

                    end

                end

                drms = rms(data(:));
                        
                allWfs(evt, ks).Z        = allWfs(evt, ks).Z/drms;
                allWfs(evt, ks).north    = allWfs(evt, ks).north/drms;
                allWfs(evt, ks).east     = allWfs(evt, ks).east/drms;
           
                dN = allWfs(evt,ks).north;
                dE = allWfs(evt,ks).east;

                allWfs(evt,ks).R = cosd(Az)*dN + sind(Az)*dE;
                allWfs(evt,ks).T = -sind(Az)*dN + cosd(Az)*dE;

                r                   = abs(hilbert(sqrt(allWfs(evt, ks).north.^2 + allWfs(evt, ks).east.^2)));
                ind                 = find(t > pre(evt), 1);

                allWfs(evt, ks).snr = max(r)/(2*std(r(1:ind)) + mean(r(1:ind)));

            end

            myTrace = [];
           

        end
       
    end

    %remap for names that need to be changed
    for ks = 1:length(S)

        check = strcmp(S(ks).StationCode, stations_to_change);

        if ~any(check)

            continue

        end

        ksnew = find(strcmp(statiion_new_name(check), {S.StationCode}));%an index

        for k = 1:length(E)

            if any(isnan(allWfs(k, ksnew).station)) && any(~isnan(allWfs(k, ks).station))%any

                tmp              = allWfs(k, ksnew);
                allWfs(k, ksnew) = allWfs(k, ks);
                allWfs(k, ks)    = tmp;

            end

        end

    end

    %remove station that don't have anythihg for
    for ks = 1:length(S)

        kill(ks) = 0;

        %keep it if you have anything for any event
        if any(~isnan([allWfs(:, ks).latitude]))

            continue

        end

        kill(ks) = 1;

    end

    if any(kill)

        allWfs(:,logical(kill)) = [];

    end

    %remove events that don't have anythihg for
    for es = 1:length(E)

        kill(es) = 0;

        %keep it if you have anything for any event
        if any(~isnan([allWfs(es, :).latitude]))

            continue

        end

        kill(es) = 1;

    end

    if any(kill)

        allWfs(logical(kill), :) = [];

    end

end