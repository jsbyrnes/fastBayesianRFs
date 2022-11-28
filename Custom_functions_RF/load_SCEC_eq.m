function allWfs = load_SCEC_eq(Parameters, S, startTime, window, lat_evt, lon_evt)

    %some networks are right handed, but typical convections are left
    %handed. So for some listed here you have to flip stuff to make it work

    left_handed_networks = { 'YL' };

    stations_to_change = { 'G02B', 'G04B'};
    statiion_new_name  = { 'G02',  'G04' };

    %first, how much padding to do to account for the 50% taper. 

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
        
    t = Parameters.t;    

    for ks = 1:length(S)

        %first, scale the array
        allWfs(1, ks).north         = nan;
        allWfs(1, ks).east          = nan;
        allWfs(1, ks).Z             = nan;
        allWfs(1, ks).R             = nan;
        allWfs(1, ks).T             = nan;
        allWfs(1, ks).latitude      = nan;
        allWfs(1, ks).longitude     = nan;
        allWfs(1, ks).snr           = nan;
        allWfs(1, ks).station       = nan;
        allWfs(1, ks).phase_list    = nan;
        allWfs(1, ks).phase_times   = nan;
        allWfs(1, ks).delta         = nan;
        allWfs(1, ks).backazi       = nan;
        allWfs(1, ks).orientation   = nan;
        allWfs(1, ks).orientation_k = nan;
        allWfs(1,ks).t0             = nan;
        
        %first check to see if this station was active when the earthquake
        %happened
        
        startTime = datenum(startTime);        
        endTime   = datenum(startTime) + window/(60*60*24);

        startTime = datestr(startTime,'yyyy-mm-dd HH:MM:SS.FFF');   
        endTime   = datestr(endTime,'yyyy-mm-dd HH:MM:SS.FFF');
    
        try
                        
            myTrace = irisFetch.Traces(S(ks).NetworkCode,...
                S(ks).StationCode,'*',channel_string,startTime, endTime);
                                    
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

        if ~isempty(myTrace)
                                                   
            %upper limits    
            myTrace = wfButterworth(myTrace, [ Parameters.high_pass Parameters.low_pass ]);
            myTrace = wfResample_jsb(myTrace, Parameters);

            names = {myTrace.channel};
            
            %find the north channel (might be prerotation)
            vind = find(contains(names, Parameters.vertical_names),1);

            %find the north channel (might be prerotation)
            nind = find(contains(names, Parameters.north_names),1);

            %find the east channel (might be prerotation)
            eind = find(contains(names, Parameters.east_names),1);
            
            if ~isempty(vind) && ~isempty(nind) && ~isempty(eind)

                %find the one that matches
                %put in the right index                   

                allWfs(1, ks).Z     = myTrace(vind).data;

                if contains(myTrace(1).network, left_handed_networks)

                    allWfs(1, ks).north     = myTrace(eind).data;
                    allWfs(1, ks).east      = myTrace(nind).data;

                else

                    allWfs(1, ks).north     = myTrace(nind).data;
                    allWfs(1, ks).east      = myTrace(eind).data;

                end

                allWfs(1, ks).latitude  = myTrace(eind).latitude;
                allWfs(1, ks).longitude = myTrace(eind).longitude;
                allWfs(1, ks).station   = myTrace(eind).station;

                [D, Az] = distance(allWfs(1, ks).latitude, allWfs(1, ks).longitude, lat_evt, lon_evt);

                allWfs(1, ks).delta       = D;%in secs into trace
                allWfs(1, ks).backazi     = Az;%in secs into trace
                
                data = [myTrace(vind).data myTrace(nind).data myTrace(eind).data];%for normalization

                if wrapTo360(myTrace(eind).azimuth - myTrace(nind).azimuth) == 90

                    dN = allWfs(1,ks).north;
                    dE = allWfs(1,ks).east;

                    allWfs(1,ks).north = cosd(myTrace(nind).azimuth)*dN - sind(myTrace(nind).azimuth)*dE;
                    allWfs(1,ks).east  = sind(myTrace(nind).azimuth)*dN + cosd(myTrace(nind).azimuth)*dE;

                end

            else

                continue

            end

            if wrapTo360(myTrace(eind).azimuth - myTrace(nind).azimuth) == 90

                allWfs(1,ks).orientation   = 0;%orientated station
                allWfs(1,ks).orientation_k = 1/(Parameters.orientation_std)^2;%...maybe

            else %ocean bottom seismometers are both listed as zero

                allWfs(1,ks).orientation   = 0;%number doesn't matter
                allWfs(1,ks).orientation_k = 1/(100)^2;%basically a uniform distribution

            end

            drms = rms(data(:));
                    
            allWfs(1, ks).Z        = allWfs(1, ks).Z/drms;
            allWfs(1, ks).north    = allWfs(1, ks).north/drms;
            allWfs(1, ks).east     = allWfs(1, ks).east/drms;
       
            dN = allWfs(1,ks).north;
            dE = allWfs(1,ks).east;

            %R = -1*(cosd(H.BAZ)*myTrace(nind).data + sind(H.BAZ)*myTrace(eind).data);%%%%???????
            %T = -1*(-sind(H.BAZ)*myTrace(nind).data + cosd(H.BAZ)*myTrace(eind).data);

            allWfs(1,ks).R = -1*(cosd(Az)*dN + sind(Az)*dE);
            allWfs(1,ks).T = -1*(-sind(Az)*dN + cosd(Az)*dE);

            %get some arrival times
            [t0,t1,y0,y1,dt0,E,h_E] = envelopePick(t,allWfs(1,ks).Z,allWfs(1,ks).delta/0.481,5,5,1,10,100,1e-3);
            allWfs(1,ks).t0 = t0;
            
        end

        myTrace = [];
       

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

end