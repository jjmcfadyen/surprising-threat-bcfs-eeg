function sim_rt = vis_lba(nTrials,params,rt,rt_lims)

sim_plot = false;

tau = .001; % time resolution

v = params(1);
A = params(2);
b = params(3);
sv = params(4);
t0 = params(5);

min_rt = rt_lims(1);
max_rt = rt_lims(2);

if sim_plot
    figure
    xidx = 0:tau:max_rt;
end

sim_rt = [];
progress_idx = round(linspace(1,nTrials,20));
for trl = 1:nTrials
    
    if trl == progress_idx(1,1)
        disp(['Simulating data... ' num2str(round(trl/nTrials)*100,2) '%'])
        progress_idx(:,1) = [];
    end
    
    % Get starting point (somewhere between 0 and A)
    this_start = rand.*A;
    
    % Get drift rate
    this_v = normrnd(v, sv);
    
    % Get time to threshold
    t = (b-this_start)./this_v;
    
    % Add on non-decision time
    r = t0 + t;
    
    if r > max_rt; r = max_rt; end
    if r < min_rt; r = 0; end
    
    % Plot
    if sim_plot
        scatter(0,this_start,30,'filled','k'); hold on % Starting point
        P = plot([0 length(0:tau:t0) length(0:tau:r)]*tau, ...
             [this_start this_start b],'k');
        P.Color(4) = .1; % alpha
        plot([0 max_rt],[b b],'r') % response threshold
        scatter(r,b,30,'filled','b') % RT
        xlabel('Time (seconds)')
        ylabel('Evidence Accumulated')
        xlim([0 max_rt])
    end
    
    sim_rt(trl) = r;
     
end
if sim_plot
    plot([mean(sim_rt) mean(sim_rt)],[0 b],'LineWidth',2,'Color','b','LineStyle','--')
end
disp(['Mean simulated RT = ' num2str(mean(sim_rt)) ]);

figure; histogram(sim_rt); title('Simulated Data'); xlabel('RT (seconds)'); hold on
ax = gca;
plot([mean(sim_rt) mean(sim_rt)],[ax.YLim],'Color','b','LineWidth',2,'LineStyle','--')
plot([mean(rt) mean(rt)],[ax.YLim],'Color','r','LineWidth',2,'LineStyle','--')
legend({'Simulated RTs','Mean Simulated RT','Mean Observed RT'})