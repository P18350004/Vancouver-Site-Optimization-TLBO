clc
clear all
tic
%%
lb = [0 0 0 0 0 0 0 0 0 70 70 70 70 70 70 70 70 70];              %lower bound
ub = [1 2 3 4 5 6 8 10 14 200 200 200 200 200 200 200 200 200];       %upper bound
prob = @Fitness_misfit;                          %fitness function

%%
Np = 200;   % population size
T = 100;   % number pf iterations
%%

f = NaN(Np,1); %% to store fitness function value

BestFitIter = NaN(T+1,1);   

D = length(lb);   % to determine the no of decision variables

P = repmat(lb,Np,1) + repmat((ub-lb),Np,1).*rand(Np,D);  % generation of initial population

for p =1:Np
    f(p) = prob(P(p,:));     % evaluating the fitness function of initial population
end

BestFitIter(1) = min(f);

%%
for t = 1:T
     
    for i = 1:Np
        %% teacher phase
        Xmean = mean(P);  % determine the mean of population
        
        [~,ind] = min(f);  % finding the loacation of teacher
        
        Xbest = P(ind,:);   % copying the solution acting as teacher
        
        TF = randi([1 2],1,1);  % Teaching factor
        
        Xnew = P(i,:) + rand(1,D).*(Xbest - TF*Xmean);
        
        Xnew = min(ub,Xnew);
        Xnew = max(lb, Xnew);
        
        fnew = prob(Xnew);
        
        if (fnew < f(i))
            P(i,:) = Xnew;
            f(i) = fnew;
        end
      %% learner phase  
        p = randi([1 Np],1,1);
        
        

        while i == p
            p = randi([1 Np],1,1);
        end
        
        if f(i) < f(p)
            
            Xnew = P(i,:) + rand(1,D).*(P(i,:) - P(p,:));
        else
            Xnew = P(i,:) - rand(1,D).*(P(i,:) - P(p,:));
        end
        
        Xnew = min(ub,Xnew);
        Xnew = max(lb, Xnew);
        
        fnew = prob(Xnew);
        
        if (fnew < f(i))
            P(i,:) = Xnew;
            f(i) = fnew;
        end
    end
    
    BestFitIter(t+1) = min(f);
    disp(['Iteration' num2str(t) ': Best Fitness =' num2str(BestFitIter(t+1))])
end

[bestfitness, ind] = min(f)
bestsol = P(ind,:);
%%
% subplot(1,2,1)
FigWidth = 12; % cm
FigHeight = 12; % cm
FigFontSize = 14; % pt

figure
plot(0:T, BestFitIter','LineWidth',3);
xlabel('Iterations','Fontname','Times New Roman');
ylabel('Best Fitness Value','Fontname','Times New Roman');
title('Teaching-learning-based Optimization','Fontname','Times New Roman');
set(gca,'Fontsize',14,'Fontname','Times New Roman');
set(gcf,'units','centimeters')
pos = [2, 2, FigWidth, FigHeight]; 
set(gcf,'Position',pos)

toc
clock
save Result
