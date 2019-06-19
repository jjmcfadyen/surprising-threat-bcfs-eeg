function [chiSquare, params] = fitModel(rt,a,params,nIterations,tau,max_rt,make_figures)

% a is fixed (boundary)

% Model parameters
Ter = params(1); % non-decision time
eta = params(2); % SD of drift rate
st = params(3); % SD of Ter
v = params(4); % mean drift rate
s = params(5); % noise in random walk

sim_rt = vis_sim(rt,[a params],nIterations,tau,max_rt,make_figures);

% Fit model
Q = .05:.05:.95;
Q_rt = quantile(rt,Q);
Q_sim = quantile(sim_rt,Q);

E = []; % Expected Values (no. of OBSERVED responses lying between the quantiles of the original RT distribution)
O = []; % Observed Values (no. of PREDICTED responses lying between the quantiles of the original RT distribution)
chiSquare = 0;
for q = 2:length(Q)
    if q == 2
        E(q-1) = (sum(sim_rt < Q_rt(q))/size(sim_rt,1))*length(rt);
        O(q-1) = (sum(rt < Q_rt(q))/size(rt,1))*length(rt);
    elseif q == length(Q)
        E(q-1) = (sum(sim_rt >= Q_rt(q-1))/size(sim_rt,1))*length(rt);
        O(q-1) = (sum(rt >= Q_rt(q-1))/size(rt,1))*length(rt);
    else
        E(q-1) = (sum(sim_rt < Q_rt(q) & sim_rt >= Q_rt(q-1))/size(sim_rt,1))*length(rt);
        O(q-1) = (sum(rt < Q_rt(q) & rt >= Q_rt(q-1))/size(rt,1))*length(rt);
    end
    
    if(E(q-1) == 0)
        chiSquare = realmax;
    else
        chiSquare = chiSquare + ((O(q-1)-E(q-1))^2)/E(q-1);
    end
end

if make_figures
    figure
    bar([E; O])
    title(['Expected vs. Observed Values (chi square = ' num2str(round(chiSquare,2)) ')'])
    set(gca,'XTickLabels',{'Expected','Observed'})
    ylabel('Count')

    disp(['**Chi Square = ' num2str(round(chiSquare,2))])
end

end