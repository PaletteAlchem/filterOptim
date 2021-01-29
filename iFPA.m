function [fitness_history, best_solution, best_fitness] = iFPA(obj_function, n, Lb, Ub, max_iteration, p, add_initial_sol)
%% Algorithm Parameters;

dim = length(Lb);
t = 1; % iteration index;
%p = 0.8; % p is the change probability;

%% Initialization;
% allocate memory;
Sol(n,dim) = 0;
Fitness(n) = 0;
for i = 1:n
    % generate n particles and evaluate their fitness;
    Sol(i,:) = Lb + (Ub-Lb).*rand(1,dim);
    Fitness(i) = obj_function(Sol(i,:));
end
Sol(n+1,:) = add_initial_sol;
Fitness(i+1) = obj_function(Sol(n+1,:));

% Find global best solution and fitness;
[fmin, Imin] = min(Fitness);
best_solution = Sol(Imin, :);
best_fitness = fmin;
%phi1 = 1;
fitness_history(max_iteration) = 0;
fitness_history(1) = min(Fitness(1:end-1));
%% Mainloop
while t <= max_iteration
    for i = 1:(n+1)
        if rand > p
            L = levy(dim);
            dS = L.*(Sol(i,:)-best_solution);
            S = Sol(i,:) + dS;
            if any(S>=Ub) || any(S<=Lb)
                newS = Lb + (Ub-Lb).*rand(1,dim);
                S(S>=Ub) = newS(S>=Ub);
                S(S<=Lb) = newS(S<=Lb);
            end
        
        else
            U = rand;
            J = rand;
            JK = randperm(n);
            S = Sol(i,:) + U * (Sol(JK(1),:)-Sol(JK(2),:))+ J * (best_solution-Sol(JK(3),:));
            if any(S>=Ub) || any(S<=Lb)
                newS = Lb + (Ub-Lb).*rand(1,dim);
                S(S>=Ub) = newS(S>=Ub);
                S(S<=Lb) = newS(S<=Lb);
                
            end
        end
        %S = checkbounds(S, Lb, Ub);
        
        
        
        Fnew = obj_function(S);
        if Fnew < Fitness(i)
            Fitness(i) = Fnew;
            Sol(i,:) = S;
            if Fnew < best_fitness
                best_fitness = Fnew;
                best_solution = S;
            end
        end
    end
    fitness_history(t+1) = best_fitness;
    t = t + 1;
    disp(['At ',num2str(t), 'th iteration, best solution:', num2str(best_fitness)]);
end
end


function L = levy(d)
%%
gamma = 0.1;
beta = 1.5;
sigma = 0.7;
u = randn(1,d)*sigma;
v = randn(1,d);
step = u./abs(v).^(1/beta);
L = gamma.*step;
end
