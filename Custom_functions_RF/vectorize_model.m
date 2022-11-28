function model = vectorize_model(model, Parameters)
    %non-dimensionalized model vector 

    model.vector = [];

    [nevt, nsta] = size(model.sig);
    dt      = model.t(2) - model.t(1);
    minf    = 1/(model.t(end)*0.5);
    maxf    = 0.5/dt;

    for k = 1:model.n

        model.vector = [ model.vector; model.A{k}];

    end

    model.vector = [ model.vector; model.t0];
    model.vector = [ model.vector; model.w];

    %error terms
    model.vector = [ model.vector; (model.sig(:) - Parameters.sig_range(1))/Parameters.sig_range(2) ];
    %model.vector = [ model.vector; (model.sigT(:) - Parameters.sig_range(1))/Parameters.sig_range(2) ];
     
    if Parameters.use_covarience

        model.vector = [ model.vector; (model.r - Parameters.r_range(1))/Parameters.r_range(2) ];
        model.vector = [ model.vector; (model.f - Parameters.f_range(1))/Parameters.f_range(2) ];
        %model.vector = [ model.vector; model.a];

    end

    model.nparam = length(model.vector);

end

%     %source parameters
% 
%     if nsta > 1
% 
%         for k = 1:nevt
%     
%             model.vector = [ model.vector; model.sources{k} ];
%     %         model.vector = [ model.vector; model.sA{k} ];
%     %         model.vector = [ model.vector; model.sB{k} ];
%     %         model.vector = [ model.vector; model.st{k}/model.t(end) ];
%     %         model.vector = [ model.vector; (model.sf{k} - minf)/(maxf - minf) ];
%     
%         end
% 
%     end
