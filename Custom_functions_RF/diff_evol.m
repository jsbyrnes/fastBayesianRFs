function BestSol = diff_evol(CostFunction,nPop,lb,ub, TolX)

%%%%%%
%jsb rewrite of 
%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPEA107
% Project Title: Implementation of Differential Evolution (DE) in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%

    VarSize = length(lb);

    beta_min = 0.2;   % Lower Bound of Scaling Factor
    beta_max = 0.8;   % Upper Bound of Scaling Factor
    
    pCR = 0.2;        % Crossover Probability

    empty_individual.Position = [];
    empty_individual.Cost     = [];
    
    BestSol.Cost = inf;
    
    pop = repmat(empty_individual,nPop,1);
    
    for i=1:nPop
    
        pop(i).Position = rand(VarSize, 1).*(ub - lb) + lb;
        pop(i).Cost     = CostFunction(pop(i).Position);
        
        if pop(i).Cost<BestSol.Cost

            BestSol=pop(i);
        
        end
        
    end
    
    %% DE Main Loop
    
    for it = 1:1e9
        
        for i = 1:nPop
            
            x = pop(i).Position;
            A = randperm(nPop);
            A(A==i)=[];
            
            a = A(1);
            b = A(2);
            c = A(3);
            
            % Mutation
            %beta=unifrnd(beta_min,beta_max);
            %beta = unifrnd(beta_min,beta_max,VarSize);
            beta = rand(1, VarSize)*(beta_max - beta_min) + beta_min;
            y    = pop(a).Position+beta.*(pop(b).Position-pop(c).Position);
            y    = max(y, lb);
		    y    = min(y, ub);
		    
            % Crossover
            z  = zeros(size(x));
            j0 = randi([1 numel(x)]);
            
            for j = 1:numel(x)
                
                if j==j0 || rand<=pCR
                
                    z(j)=y(j);
                
                else
                
                    z(j)=x(j);
                
                end
            
            end
            
            NewSol.Position = z;
            NewSol.Cost     = CostFunction(NewSol.Position);
            
            if NewSol.Cost<pop(i).Cost
            
                pop(i)=NewSol;
                
                if pop(i).Cost<BestSol.Cost
                
                    BestSol=pop(i);
                
                end
            
            end
            
        end
                    
        % Show Iteration Information
        disp(['Iteration ' num2str(it) ': Best Cost,' num2str(BestSol.Cost) ...
            ' with tolerence of ' num2str(range([pop(:).Cost]))]);

        if range([pop(:).Cost]) < TolX

            break

        end

    end
    
end