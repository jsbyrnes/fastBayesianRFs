function [irF, irS] = ir_rotating_layer(delay, phi, dphi, tSA, om, polarization, prelogf)
    %phi is the mean angle. dphi is deviation at top and bottom. Angles are
    %positive clockwise from north, with positive deviation being positive
    %rotation towards the surface (so rotating from 20 at the surface to
    %30 at the base is phi=25 and dphi = -5, rotating from 30 to 20 is dphi = 5);
    %coordiantes gives the orientation of component 1 (2 is +90)

    %impulse responses are in thier polarization direction (irF is in the
    %fast direction and irS is in the slow direction). Rotated elsewhere;
    %this is necessary to calculate the impulse response in the second
    %layer.

    fft_scale = length(om)*(om(2) - om(1))/(2*pi);%just the sampling frequency

    %first, get the single - scattered waves
    %fast
    scale = cos(phi - dphi - polarization);
    irF   = real(ifft(scale*exp(1i*om*delay/2)));%.*[cosd(thetahat(1) + 90*1^key(end)); sind(thetahat(end) + 90*1^key(1)) ]));

    %slow
    scale = cos(phi - dphi - polarization + pi/2);
    irS   = real(ifft(scale*exp(-1i*om*delay/2)));%.*[cosd(thetahat(1) + 90*1^key(end)); sind(thetahat(end) + 90*1^key(1)) ]));    
    
    if abs(dphi) > 0 && delay > 0

        idealN = delay*fft_scale;%space out deltas to the sampling rate
    
        dp = 2*dphi/idealN;
        %box car function is sinc in fourier domain. 
        phase = fft_scale*real(ifft(delay*sinc(delay*om/(2*pi))));

        %start slow and scatter to fast....
        scale = dp*cos(phi - dphi - polarization + pi/2);
        irF = irF + scale*phase;

        %now the fast to slow
        scale = -1*dp*cos(phi - dphi - polarization);
        irS = irS + scale*phase;

    end

    %attentute one of the two components
    if abs(tSA) > 0

        f = om/(2*pi);

        Amp = exp(-f*abs(tSA)/2);                 % amplitude spectrum
        %%%%need to normalize amp
        Dt  = -(abs(tSA)/pi)*prelogf; % dummy F=0, OK on next line
        Dt  = Dt - min(Dt);
        Dph = f.*(Dt);                   % delta phase shift spectrum

        if tSA > 0 %apply to the slow component

            irS = real(ifft(fft(irS).*Amp.*exp(1i*Dph)));

        elseif tSA < 0 %apply (absolute value) to the fast component

            irF = real(ifft(fft(irF).*Amp.*exp(1i*Dph)));

        end

    end

end

% 
%     if abs(dphi) > 0
% 
%         idealN = delay*fft_scale;%space out deltas to the sampling rate
% 
%         %in test, 10 layers always seems fine - even for very large
%         %rotations
% %         n = 1000;
% %  
% %         dp = 2*dphi*pi/180/n;
% % 
% %         %question - can this for loops be replaces with a box car summation?
% % 
% %         %doubly scattered from slow to fast (once at the bottom, and then
% %         %internally). Note that scale is invarient. 
% %         for s = 1:n
% % 
% %             dtsum = 0.5*(delay/n)*(s - (n - s));
% %             %do all of the projections
% %             scale = dp*cosd(phi - dphi - polarization + 90);
% %               
% %             irF = irF + real(ifft(scale*exp(-1i*om*dtsum)));
% % 
% % %             phase = real(ifft(scale*exp(-1i*om*dtsum)));%.*[cosd(thetahat(1) + 90*1^key(end)); sind(thetahat(end) + 90*1^key(1)) ]));
% % %     
% % %             irF(1, :) = irF(1, :) + phase*cosd(coordinates);
% % %             irF(2, :) = irF(2, :) + phase*sind(coordinates);
% % 
% %         end
%     
%         dp = 2*dphi*pi/180/idealN;
%         %box car function is sinc in fourier domain. 
%         phase = fft_scale*real(ifft(delay*sinc(delay*om/(2*pi))));
%         scale = dp*cosd(phi - dphi - polarization + 90);
%         irF = irF + scale*phase;
% 
%         %doubly scattered from slow to fast (once at the bottom, and then
%         %internally). Note that negative sign. 
% 
% %         n = 1000;
% %         dp = 2*dphi*pi/180/n;        
% %         
% %         for s = 1:n
% %     
% %             dtsum = 0.5*(delay/n)*(-s + (n - s));
% %             %do all of the projections
% %             scale = -1*dp*cosd(phi - dphi - polarization);
% %             
% %             irS = irS + real(ifft(scale*exp(-1i*om*dtsum)));
% % 
% % %             phase = real(ifft(scale*exp(-1i*om*dtsum)));%.*[cosd(thetahat(1) + 90*1^key(end)); sind(thetahat(end) + 90*1^key(1)) ]));
% % %     
% % %             irS(1, :) = irS(1, :) + phase*cosd(coordinates + 90);
% % %             irS(2, :) = irS(2, :) + phase*sind(coordinates + 90);
% %     
% %         end
% 
%         %dp = 2*dphi*pi/180/idealN;
%         %box car function is sinc in fourier domain. 
%         %phase = fft_scale*real(ifft(delay*sinc(delay*om/(2*pi))));
% 
%         scale = -1*dp*cosd(phi - dphi - polarization);
%         irS = irS + scale*phase;
% 
%     end
% 

        %in test, 10 layers always seems fine - even for very large
        %rotations
%         n = 1000;
%  
%         dp = 2*dphi*pi/180/n;
% 
%         %question - can this for loops be replaces with a box car summation?
% 
%         %doubly scattered from slow to fast (once at the bottom, and then
%         %internally). Note that scale is invarient. 
%         for s = 1:n
% 
%             dtsum = 0.5*(delay/n)*(s - (n - s));
%             %do all of the projections
%             scale = dp*cosd(phi - dphi - polarization + 90);
%               
%             irF = irF + real(ifft(scale*exp(-1i*om*dtsum)));
% 
% %             phase = real(ifft(scale*exp(-1i*om*dtsum)));%.*[cosd(thetahat(1) + 90*1^key(end)); sind(thetahat(end) + 90*1^key(1)) ]));
% %     
% %             irF(1, :) = irF(1, :) + phase*cosd(coordinates);
% %             irF(2, :) = irF(2, :) + phase*sind(coordinates);
% 
%         end

        %doubly scattered from slow to fast (once at the bottom, and then
        %internally). Note that negative sign. 

%         n = 1000;
%         dp = 2*dphi*pi/180/n;        
%         
%         for s = 1:n
%     
%             dtsum = 0.5*(delay/n)*(-s + (n - s));
%             %do all of the projections
%             scale = -1*dp*cosd(phi - dphi - polarization);
%             
%             irS = irS + real(ifft(scale*exp(-1i*om*dtsum)));
% 
% %             phase = real(ifft(scale*exp(-1i*om*dtsum)));%.*[cosd(thetahat(1) + 90*1^key(end)); sind(thetahat(end) + 90*1^key(1)) ]));
% %     
% %             irS(1, :) = irS(1, :) + phase*cosd(coordinates + 90);
% %             irS(2, :) = irS(2, :) + phase*sind(coordinates + 90);
%     
%         end

