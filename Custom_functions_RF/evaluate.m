function model = evaluate(model, allWfs, Parameters)
    
    [~, nsta] = size(allWfs);

    model = apply_model(model, allWfs);

    model.llh = -0.5*log(2*pi)*length(model.t) - ...
        (length(allWfs)*model.logdet + 2*length(model.t)*sum(model.sig(:)))/2 ...
       - model.phi/2;
   
    model.lp = -0.5*sum( ((model.sig(:) - Parameters.sig_range(1)).^2)/(Parameters.sig_range(2)^2));
    %model.lp = sum( ((model.sigT(:) - Parameters.sig_range(1)).^2)/(Parameters.sig_range(2)^2));
                     
    %enforce bounds on time, etc

    n = 0;
    for k = 1:length(model.t0)

         model.lp = model.lp - sum(1e9*(model.t0<0));
%         model.lp = model.lp - sum(1e9*(model.f0{k}<minf));
%         model.lp = model.lp - sum(1e9*(model.f0{k}>maxf));
%        model.lp = model.lp - sum(exp((model.t0 - model.t(end))));

        %t0       = model.t0;
        model.lp = model.lp - sum(model.t0(model.t0>0)/(model.t(end)*0.25));

        model.lp = model.lp - sum(1e9*(model.w<(2*log(model.t(2)))));
        model.lp = model.lp - sum(1e9*(model.w>(log(model.t(end)/20))));
        %model.lp = model.lp - sum(model.w.^2)/(2*0.5^2);
        %model.lp = model.lp - sum(model.w.^2)/(2*0.5^2);

        for i = 1:length(model.A{k})

            model.lp = model.lp - sum(model.A{k}.^2)/(2*0.1^2);

            n = n + length(model.A{k});

        end

    end

    if Parameters.use_covarience

        model.lp = model.lp - (((model.r - Parameters.r_range(1))^2)/Parameters.r_range(2)^2 + ...
            sum(((model.f - Parameters.f_range(1)).^2)./Parameters.f_range(2)^2))/2;
        %model.lp = model.lp - 1e9*sum(model.a<0);
        %model.lp = model.lp - 1e9*sum(model.a>1);

        if exp(model.f) > Parameters.low_pass
            
            model.lp = -1e9;
    
        end

    end

    model.lpst = model.llh + model.lp;
    model.BIC  = -2*model.lpst + log(length(model.R)*length(model.t))*(2*sum(model.n) + n);

end

%%%%%%%%
%old search for alignment
%                     lags = 1:length(dN);
%                     for j = 1:length(lags)
% 
%                         %comutationally intensive
%                         phi(j) = md(circshift(No, lags(j)), circshift(Eo, lags(j)));
%     
%                     end
%     
%                     [R(2), idx] = min(phi);
%                     %find neighbors, assuming periodic
%                     if idx==1
%     
%                         R(1)=phi(end);
%                         R(3)=phi(2);
%     
%                     elseif idx==numel(phi)
%     
%                         R(1)=phi(end-1);
%                         R(3)=phi(1);
%     
%                     else
%     
%                         R(1)=phi(idx-1);
%                         R(3)=phi(idx+1);
%     
%                     end
%     
%                     c     = (R(3)-R(1))/(2*(2*R(2)-R(1)-R(3)));
%                     %lag   = mod(idx-1+floor(length(lags)/2),length(lags))-floor(length(lags)/2);%integer part
%                     delay = idx+c;%delay estimate
%                     N     = delay_continuous(No, Parameters.sample_rate, delay/Parameters.sample_rate );
%                     E     = delay_continuous(Eo, Parameters.sample_rate, delay/Parameters.sample_rate );


%double grid search for amplitude
%                     phi = []
%                     amps = (-6:0.05:3)
%                     for j = 1:length(amps)
%     
%                         phi(j) = md(exp(amps(j))*N, exp(amps(j))*E);
%     
%                     end
%     
%                     [mphi, idx] = min(phi);
% 
%                     if idx > 1
% 
%                         %refined search - needs to be extremely accurate
%                         phi = [];
%                         amps = (-0.5:0.05:0.5) + amps(idx);
%                         for j = 1:length(amps)
%         
%                             %computationally intensive
%                             phi(j) = md(exp(amps(j))*N, exp(amps(j))*E);
%         
%                         end
%                         [mphi, idx] = min(phi);        
%                         c     = (phi(idx+1) - phi(idx-1))/(2*(2*mphi - phi(idx-1) - phi(idx+1)));
%                         %ampi  = mod(idx - 1 + floor(length(amps)/2), length(amps)) - floor(length(amps)/2);%integer part
%                         amp   = amps(idx) + c*0.05;%amp estimate
% 
%                     else
% 
%                         amp = -6;%basically zero, its a huge mismatch. 
% 
%                     end
% 
%                     N = exp(amp)*N;
%                     E = exp(amp)*E;

%         if sum(model.dt(k, :)) > model.dt(k, 1)%check for a second layer
% 
%             %account for rotation
% %             dphi     = abs( (model.fast_dir(k, 1) + model.fast_dir_rotation(k,1)) ...
% %                 - (model.fast_dir(k, 2) - model.fast_dir_rotation(k,2)));
%             dphi     = abs( model.fast_dir(k, 1) - model.fast_dir(k, 2));
% 
%             if dphi > pi/2
% 
%                 dphi = dphi - pi/2;
% 
%             end
% 
%             dphi = dphi/(pi/2);
% 
%             model.lp = model.lp + (Parameters.beta - 1)*(log(dphi) + log(1 - dphi));
% 
%         end

%     if nsta > 1
% 
%         for k = 1:length(model.sources)
%     
%             %model.lp = model.lp - 0.5*sum(model.sources{k}.^2);
%     
%     %         model.lp = model.lp - 0.5*sum(model.sA{k}.^2);
%     %         model.lp = model.lp - 0.5*sum(model.sB{k}.^2);
%     % 
%     %         model.lp = model.lp - sum(exp(-model.st{k}));
%     %         model.lp = model.lp - sum(exp((model.st{k} - model.t(end))));
%     %         model.lp = model.lp - sum(exp(-(model.sf{k} - minf)));
%     %         model.lp = model.lp - sum(exp((model.sf{k} - maxf)));
%     
%             ns = ns + length(model.sources{k});
%     
%         end
% 
%     end

