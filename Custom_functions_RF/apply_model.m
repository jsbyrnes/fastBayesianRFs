function model = apply_model(model, allWfs)
%the last four optional outputs are really only useful when evt_ind and
%sta_ind are scalars

    [nevt, nsta] = size(allWfs);
    dt = model.t(2);

    %make the receiver function

    if model.update_rf

        rf = zeros(size(model.t));
            
        for i = 1:model.n
            
            w = exp(model.w(i));
    
            pulse = exp(-((model.t - model.t0(i)).^2)/(2*w^2));
    
            A = model.A{i};
    
            for j = 1:length(A)
    
                rf    = rf + A(j)*pulse;
                pulse = gradient(pulse, dt);
                pulse = pulse/max(abs(pulse));
    
            end
            
        end
    
        model.rf = rf;
    
        rf = fft(rf);
    
        %convolve
        for k = 1:nevt
                            
            source = fft(allWfs(k).Z);
            R      = real(ifft(rf.*source));
            
            %model.Z{k,kk} = Z;
            model.R{k} = R;
            %model.T{k,kk} = T;
            
        end

    end

    for k = 1:nevt

        model.phi_vec(k) = (model.R{k} - allWfs(k).R)'*(model.Cinv/...
            (exp(model.sig(k))^2))*(model.R{k} - allWfs(k).R);
    end

    model.phi     = sum(model.phi_vec(:));

end


%     dt = model.t(2) - model.t(1);
% 
%     [nevt, nsta] = size(allWfs);
% 
%     for k = evt_ind
% 
%         %now build the source
% 
%         if nsta > 1
% 
%             source = harmonic_sum(model.sources{evt_ind}, length(model.t));
%             %source = sum(cell2mat(arrayfun( @(A,B,t,f) returngabor(A,B,t,f, model.t), ...
%             %    model.sA{evt_ind}, model.sA{evt_ind}, model.st{evt_ind}, model.sf{evt_ind}, 'UniformOutput', false)')')';
%             source = fft(source);
% 
%         end
% 
%         for kk = sta_ind
% 
%             if isnan(allWfs(k, kk).north)
% 
%                 model.phi_vec(k, kk) = 0;
%                 model.sig(k, kk)     = 0;
% 
%             else
%                    
%                 if nsta > 1
% 
%                     Z = zeros(size(model.t));
%                     R = zeros(size(model.t));
%                     T = zeros(size(model.t));
%     
%                     az = model.az{kk};
%                     ar = model.ar{kk};
%                     at = model.at{kk};
%                     t0 = model.t0{kk};
%                     w  = exp(model.w{kk});
%                     
%                     for i = 1:model.nrf(kk)
%                 
%                         pulse = exp(-((model.t - t0(i)).^2)/(2*w(i)^2));
%                         Z = Z + az(i)*pulse;
%                         R = R + ar(i)*pulse;
%                         T = T + at(i)*pulse;
%                 
%                     end
%     
%                     model.Zrf{k,kk} = Z;
%                     model.Rrf{k,kk} = R;
%                     model.Trf{k,kk} = T;
%     
%                     Z = real(ifft(fft(Z).*source));
%                     R = real(ifft(fft(R).*source));
%                     T = real(ifft(fft(T).*source));
%     
%                     model.Z{k,kk} = Z;
%                     model.R{k,kk} = R;
%                     model.T{k,kk} = T;
%     
%                     model.phi_vec(k, kk) = (Z - allWfs(k,kk).Z)'*(model.Cinv/...
%                         (exp(model.sig(k, kk))^2))*(Z - allWfs(k,kk).Z) ...
%                         + (R - allWfs(k,kk).R)'*(model.Cinv/...
%                         (exp(model.sig(k, kk))^2))*(R - allWfs(k,kk).R) ...
%                         + (T - allWfs(k,kk).T)'*(model.Cinv/(exp(model.sig(k, kk))^2))...
%                         *(T - allWfs(k,kk).T);
%     
%                 else
% 
%                     Z = zeros(size(model.t));
%                     R = zeros(size(model.t));
%                     T = zeros(size(model.t));
%     
%                     Ar = model.Ar{kk};
%                     Br = model.Br{kk};
%                     At = model.At{kk};
%                     Bt = model.Bt{kk};
% %                     A0 = model.A0{kk};
% %                     B0 = model.B0{kk};
% %                     A1 = model.A1{kk};
% %                     B1 = model.B1{kk};
% %                     A2 = model.A2{kk};
% %                     B2 = model.B2{kk};
%                     t0 = model.t0{kk};
%                     f0 = model.f0{kk};
%                     dp = model.dp{kk};
%                     w  = exp(model.w{kk});
%                     
%                     for i = 1:model.nrf(kk)
%             
% %                         aziterm =  + A1(i)*cosd(allWfs(k,kk).backazi) + B1(i)*sind(allWfs(k,kk).backazi) ...
% %                              + A2(i)*cosd(2*allWfs(k,kk).backazi) + B2(i)*sind(2*allWfs(k,kk).backazi);
% 
%                         pulse = exp(-((model.t - t0(i)).^2)/(2*w(i)^2));
%                         %R = R + returngabor(Ar(i), Br(i),t0(i),f0(i),w(i),model.t);%ar(i)*pulse;
%                         R = R + Ar(i)*pulse;
%                         T = T + At(i)*pulse;
% 
% %                         aziterm =  + A1(i)*cosd(allWfs(k,kk).backazi + 90) + B1(i)*sind(allWfs(k,kk).backazi + 90) ...
% %                              + A2(i)*cosd(2*(allWfs(k,kk).backazi + 90)) + B2(i)*sind(2*(allWfs(k,kk).backazi + 90));
%                         
%                         %T = T + returngabor(At(i), Bt(i),t0(i),f0(i),w(i),model.t);%at(i)*pulse;
%             
%                     end
%     
%                     model.Zrf{k,kk} = Z;
%                     model.Rrf{k,kk} = R;
%                     model.Trf{k,kk} = T;
%     
%                     source = fft(allWfs(k,kk).Z);
% 
%                     R = real(ifft(fft(R).*source));
%                     T = real(ifft(fft(T).*source));
%     
%                     model.Z{k,kk} = Z;
%                     model.R{k,kk} = R;
%                     model.T{k,kk} = T;
%     
%                     model.phi_vec(k, kk) = (R - allWfs(k,kk).R)'*(model.Cinv/...
%                         (exp(model.sigR(k, kk))^2))*(R - allWfs(k,kk).R) ...
%                         + (T - allWfs(k,kk).T)'*(model.Cinv/(exp(model.sigT(k, kk))^2))...
%                         *(T - allWfs(k,kk).T);
% 
%                 end
% 
%             end
% 
%         end
% 
%     end
