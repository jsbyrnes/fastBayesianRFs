function [index, field, scale, evtind, staind] = map_index(index0, model, allWfs, TD_parameters)
%type currently not in use

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get the linear index of what will be changed
    [nevt, nsta] = size(allWfs);
    nwavelet     = numel(model.wavelet);%n in wavelets
    nsplit       = numel(model.A);%number of splitting parameters. Could be clustered

    scale = 1;
    index = index0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %map the index to the structure fields
    evtind = [];
    staind = [];
    type   = [];
    field  = [];
    %change to the wavelet
    if index <= nwavelet

        [evtind, ~] = ind2sub(size(model.wavelet), index);
        staind      = 1:nsta;
        
        field = 'wavelet';
        type  = 'all';

    else

        index = index - nwavelet;

    end

    %change to polarization
    if index <= nevt && isempty(field)

        staind = 1:nsta; 
        evtind = index;%maps directly
        field  = 'polarization';
        type   = 'all';
        scale  = TD_parameters.polarization_std;

    elseif isempty(field)

        index = index - nevt;

    end

    %%%%amp
    if index <= nsta*nevt && isempty(field)

        [evtind, staind] = ind2sub(size(allWfs), index);
        
        field = 'amp';
        scale = 1;

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %%%%delay
    if index <= nsta*nevt && isempty(field)

        [evtind, staind] = ind2sub(size(allWfs), index);
        
        field = 'shift';
        scale = TD_parameters.Delaystd;

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %%%%tS
    if index <= nsta*nevt && isempty(field)

        [evtind, staind] = ind2sub(size(allWfs), index);
        
        field = 'dtS';
        scale = TD_parameters.dtS(2);

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %change to A/dt
    if index <= nsplit && isempty(field)

        if TD_parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end

        evtind      = 1:nevt;

        field = 'A';
        scale = TD_parameters.ABstd;%same search size for each
        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to B/fastdir
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        if TD_parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end

        evtind      = 1:nevt;
            
        field = 'B';
        scale = TD_parameters.ABstd;%same search size for each

        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to rotation
    if index <= nsplit && isempty(field)

        %A works for size on all splitting parameters
        [~, lind] = ind2sub(size(model.A), index);
        if TD_parameters.cluster

            staind = 1:nsta;

        else

            [staind, ~] = ind2sub(size(model.A), index);%dont need layer ind

        end
        evtind = 1:nevt;
        field  = 'fast_dir_rotation';
        scale = TD_parameters.rotation_std;

        type  = 'all';

    elseif isempty(field)

        index = index - nsplit;

    end

    %change to orientations

    if TD_parameters.orientations

        if index <= nsta && isempty(field)
    
            staind = index;%applies directly
            evtind = 1:nevt;
            
            field = 'sta_or';
            scale = min([ TD_parameters.sta_err(index), pi/2]);
            type  = 'all'; 
    
        elseif isempty(field)
    
            index = index - nsta;
    
        end

    end
    
    %change to sigma
    if index <= nsta*nevt && isempty(field)

        [evtind, staind] = ind2sub(size(allWfs), index);
        
        field = 'sig';
        scale = TD_parameters.sig_range(2);
        type  = 'error';%not really but this is where it is evaluated

    elseif isempty(field)

        index = index - nsta*nevt;

    end

    %change to r
    if index == 1 && isempty(field)

        index = 1;%scalar
        
        evtind = 1:nevt;
        staind = 1:nsta;
        scale = TD_parameters.r_range(2);
        field = 'r';
        type  = 'error';%not really but this is where it is evaluated

    elseif isempty(field)

        index = index - 1;

    end

    %change to f
    if index == 1 && isempty(field)

        evtind = 1:nevt;
        staind = 1:nsta;
        field  = 'f';
        scale  = TD_parameters.f_range(2);
        type   = 'error';%not really but this is where it is evaluated

    elseif isempty(field)

        index = index - 1;

    end

end