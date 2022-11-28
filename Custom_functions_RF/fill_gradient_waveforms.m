function model = fill_gradient_waveforms(model, allWfs, Parameters)

    del = zeros(length(model.vector),1);%overwrite it

    if Parameters.parallel

        parfor k = 1:length(model.vector)
    
            del(k, 1) = internal_loop(model, allWfs, Parameters, k);
    
        end

    else

        for k = 1:length(model.vector)
    
            del(k, 1) = internal_loop(model, allWfs, Parameters, k);
    
        end

    end

    model.del = del;

end

function [ del, dN, dE ] = internal_loop(model, allWfs, Parameters, index)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get the linear index of what will be changed
    [nevt, nsta] = size(allWfs);

    scale  = 1;

    index0 = index;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %map the index to the structure fields

    h = Parameters.h;

    if nsta > 1

        for k = 1:nevt
    
            if index <= length(model.sources{k})
    
                modelp = model;
                modeln = model;
                        
                sp        = model.sources{k};
                sp(index) = sp(index) + h;
                sn        = model.sources{k};
                sn(index) = sn(index) - h;
    
                modelp.sources{k} = sp;
                modeln.sources{k} = sn;
            
                modelp = evaluate(modelp, allWfs, Parameters, k, 1:nsta);
                modeln = evaluate(modeln, allWfs, Parameters, k, 1:nsta);
            
                del = (((modelp.lpst) - (modeln.lpst))/(2*h))*scale;
                return
    
            else
    
                index = index - length(model.sources{k});
    
            end
    
        end

    end

    for k = 1:nsta

        if nsta > 1

            %az
            if index <= model.nrf(k)
    
                modelp = model;
                modeln = model;
                        
                ap        = model.az{k};
                ap(index) = ap(index) + h;
                an        = model.az{k};
                an(index) = an(index) - h;
    
                modelp.az{k} = ap;
                modeln.az{k} = an;
            
                modelp = evaluate(modelp, allWfs, Parameters, 1:nevt, k);
                modeln = evaluate(modeln, allWfs, Parameters, 1:nevt, k);
            
                del = (((modelp.lpst) - (modeln.lpst))/(2*h))*scale;
                return
    
            else
    
                index = index - model.nrf(k);
    
            end

        end

        %ar
        if index <= model.nrf(k)

            modelp = model;
            modeln = model;
                    
            ap        = model.Ar{k};
            ap(index) = ap(index) + h;
            an        = model.Ar{k};
            an(index) = an(index) - h;

            modelp.Ar{k} = ap;
            modeln.Ar{k} = an;
        
            modelp = evaluate(modelp, allWfs, Parameters, 1:nevt, k);
            modeln = evaluate(modeln, allWfs, Parameters, 1:nevt, k);
        
            del = (((modelp.lpst) - (modeln.lpst))/(2*h))*scale;
            return

        else

            index = index - model.nrf(k);

        end

        %at
        if index <= model.nrf(k)

            modelp = model;
            modeln = model;
                    
            ap        = model.At{k};
            ap(index) = ap(index) + h;
            an        = model.At{k};
            an(index) = an(index) - h;

            modelp.At{k} = ap;
            modeln.At{k} = an;
        
            modelp = evaluate(modelp, allWfs, Parameters, 1:nevt, k);
            modeln = evaluate(modeln, allWfs, Parameters, 1:nevt, k);
        
            del = (((modelp.lpst) - (modeln.lpst))/(2*h))*scale;
            return

        else

            index = index - model.nrf(k);

        end

        %for time
        if index <= model.nrf(k)

            modelp = model;
            modeln = model;
                    
            tp        = model.t0{k};
            tp(index) = tp(index) + h;
            tn        = model.t0{k};
            tn(index) = tn(index) - h;

            modelp.t0{k} = tp;
            modeln.t0{k} = tn;
        
            modelp = evaluate(modelp, allWfs, Parameters, 1:nevt, k);
            modeln = evaluate(modeln, allWfs, Parameters, 1:nevt, k);
        
            del = (((modelp.lpst) - (modeln.lpst))/(2*h));
            return

        else

            index = index - model.nrf(k);

        end

        %for width
        if index <= model.nrf(k)

            modelp = model;
            modeln = model;
                    
            wp        = model.w{k};
            wp(index) = wp(index) + h;
            wn        = model.w{k};
            wn(index) = wn(index) - h;

            modelp.w{k} = wp;
            modeln.w{k} = wn;
        
            modelp = evaluate(modelp, allWfs, Parameters, 1:nevt, k);
            modeln = evaluate(modeln, allWfs, Parameters, 1:nevt, k);
        
            del = (((modelp.lpst) - (modeln.lpst))/(2*h));
            return

        else

            index = index - model.nrf(k);

        end

    end

    %change to sigma
    if index <= nsta*nevt

        [evtind, staind] = ind2sub(size(allWfs), index);

        modelp = model;
        modeln = model;
        
        h = Parameters.h*Parameters.sig_range(2);
    
        modelp.sigR(index) = modelp.sigR(index) + h;
        modeln.sigR(index) = modeln.sigR(index) - h;
    
        modelp = evaluate(modelp, allWfs, Parameters, evtind, staind);
        modeln = evaluate(modeln, allWfs, Parameters, evtind, staind);
    
        del = (((modelp.lpst) - (modeln.lpst))/(2*h))*Parameters.sig_range(2);
        return

    else

        index = index - nsta*nevt;

    end

    %change to sigma
    if index <= nsta*nevt

        [evtind, staind] = ind2sub(size(allWfs), index);

        modelp = model;
        modeln = model;
        
        h = Parameters.h*Parameters.sig_range(2);
    
        modelp.sigT(index) = modelp.sigT(index) + h;
        modeln.sigT(index) = modeln.sigT(index) - h;
    
        modelp = evaluate(modelp, allWfs, Parameters, evtind, staind);
        modeln = evaluate(modeln, allWfs, Parameters, evtind, staind);
    
        del = (((modelp.lpst) - (modeln.lpst))/(2*h))*Parameters.sig_range(2);
        return

    else

        index = index - nsta*nevt;

    end

    %change to r
    if index == 1
        
        modelp = model;
        modeln = model;
        
        h = Parameters.h*Parameters.r_range(2);
    
        modelp.r = modelp.r + h;
        modeln.r = modeln.r - h;
    
        modelp = build_C(modelp, model.t);
        modeln = build_C(modeln, model.t);

        modelp = evaluate(modelp, allWfs, Parameters, 1:nevt, 1:nsta);
        modeln = evaluate(modeln, allWfs, Parameters, 1:nevt, 1:nsta);
    
        del = (((modelp.lpst) - (modeln.lpst))/(2*h))*Parameters.r_range(2);
        return

    else

        index = index - 1;

    end

    %change to f
    modelp = model;
    modeln = model;
    
    h = Parameters.h*Parameters.f_range(2);

    modelp.f = modelp.f + h;
    modeln.f = modeln.f - h;

    modelp = build_C(modelp, model.t);
    modeln = build_C(modeln, model.t);

    modelp = evaluate(modelp, allWfs, Parameters, 1:nevt, 1:nsta);
    modeln = evaluate(modeln, allWfs, Parameters, 1:nevt, 1:nsta);

    del = (((modelp.lpst) - (modeln.lpst))/(2*h))*Parameters.f_range(2);
    return

end

%    if abs(del) > 100
% 
%        v = (-5:0.1:5)*scale;
% 
%         for k = 1:length(v)
%             mn(k) = model;
%             mn(k).(field)(index) = model.(field)(index) + v(k);
%             mn(k)                = apply_update(mn(k), field, Parameters.t);
%             mn(k)                = evaluate(mn(k), allWfs, Parameters, evtind, staind, type);
%             lpvec(k) = mn(k).lpst;
%         end
%         plot(v, lpvec)
%       
%    end

%     model.G = zeros(model.nparam, model.nparam);
% 
%     ndata = model.nparam - numel(allWfs) - 2;
%     ncvar = model.nparam - ndata;
% 
%     Gdata = zeros(ndata, ndata);
%     Gcvar = zeros(ncvar, ncvar);
% 
%     %Gdata matrix
%     for k = 1:((ndata^2 - ndata)/2)
% 
%         [k1,k2] = ind2sub(size(Gdata), k);
% 
%         if k2 > k1
% 
%             continue
% 
%         end
% 
%         %get derivatives
%         dn1 = dN{k1};
%         de1 = dE{k1};
%         dn2 = dN{k2};
%         de2 = dE{k2};
% 
%         for j = 1:length(dn1)
% 
%             [j1,j2] = ind2sub(size(allWfs), j);
% 
%             Gdata(k1,k2) = model.G(k1,k2) + (dn2{j}'*model.Cinv*dn1{j})/exp(model.sig(j1,j2))^2;
%             Gdata(k1,k2) = model.G(k1,k2) + (de2{j}'*model.Cinv*de1{j})/exp(model.sig(j1,j2))^2;
% 
%         end
% 
%     end
% 
%     %Gcvar matrix
%     %covarience parameters, do numerically
%     Cinv = model.Cinv;%take this out of the structure
%     
%     %r
%     modelp = model;
%     modeln = model;
%     
%     modelp.r = modelp.r + Parameters.h;
%     modeln.r = modeln.r + Parameters.h;
%     
%     [~, Cp] = build_C(modelp, Parameters.t);
%     [~, Cn] = build_C(modeln, Parameters.t);
%     
%     dCr = (Cp - Cn)/(2*Parameters.h);
%     
%     %f
%     modelp = model;
%     modeln = model;
%     
%     modelp.f = modelp.f + Parameters.h;
%     modeln.f = modeln.f + Parameters.h;
%     
%     [~, Cp] = build_C(modelp, Parameters.t);
%     [~, Cn] = build_C(modeln, Parameters.t);
%     
%     dCf = (Cp - Cn)/(2*Parameters.h);
% 
%     for k = 1:((ncvar^2 - ncvar)/2)
% 
%         [k1,k2] = ind2sub(size(Gcvar), k);
% 
%         if k2 > k1
% 
%             continue
% 
%         end
% 
%         if k1 == k2 && k1 <= (ncvar - 2)
% 
%             %sigma levels
%             %only non-zero on the diagonal
%         
%             sig = exp(model.sig(k1));
%             model.G(k1,k2) = nsta*(4/sig^2);
% 
%         end
% 
%         if k2 == (ncvar - 1) && k1 <= (ncvar - 2)
% 
%             sig = exp(model.sig(k1));
%             model.G(k1,k2) = nsta*(4/sig^2);
% 
%         end
% 
%     end

