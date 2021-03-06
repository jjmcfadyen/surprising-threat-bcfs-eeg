function analyse_behav

%% Experiments

for exp = 1:2
    
    disp('=====================================================================')
    disp(['=== EXPERIMENT ' num2str(exp) ' ===================================================='])
    disp('=====================================================================')

    output = preprocess_behav(exp,false);
    T = output.trialdata;
    
    % Exclude subjects based on missed trials
    subjects = unique(T.Subject);
    N = length(subjects);
    nmisses = nan(N,2);
    for s = 1:N
        nmisses(s,1) = sum(isnan(T.RT(T.Subject==s))); 
        nmisses(s,2) = sum(isnan(T.RT(T.Subject==s))) / sum(T.Subject==s);
    end
    
    excludeSubjects = find(nmisses(:,2) > .15);
    if any(excludeSubjects)
        disp(['Excluding subjects ' num2str(excludeSubjects')])
    end
    
    T(ismember(T.Subject,excludeSubjects),:) = [];

    for st = 1:3
        
        if st==1
            disp('--- ALL STANDARDS ---')
            idx = T.Acc & ~T.Outlier;
        elseif st==2
            disp('--- FIRST STANDARD ---')
            idx = T.Acc & ~T.Outlier & (contains(T.Type,'deviant') | contains(T.Type,'first'));
        elseif st==3
            disp('--- LAST STANDARD ---')
            idx = T.Acc & ~T.Outlier & (contains(T.Type,'deviant') | contains(T.Type,'last'));
        end

        lme = fitglme(T(idx,:), 'RT~Emotion*Expectation+(Emotion*Expectation|Subject)');
        em = emmeans(lme);
        UNEN = contrasts_wald(lme,em,[-1 1 0 0]);
        UFEF = contrasts_wald(lme,em,[0 0 -1 1]);
        UNEF = contrasts_wald(lme,em,[0 1 -1 0]);
        UFEN = contrasts_wald(lme,em,[-1 0 0 1]);

        disp(lme.anova)
        disp('[[ Standard vs Oddball - Within Emotion ]]')
        disp(['UN vs EN:    ' num2str(em.table.Estimated_Marginal_Mean([2 1])'),...
            '     (' num2str(round(diff(em.table.Estimated_Marginal_Mean([2 1])'),3)),...
            ')     p = ' num2str(UNEN.pVal)])
        disp(['UF vs EF:    ' num2str(em.table.Estimated_Marginal_Mean([4 3])'),...
            '     (' num2str(round(diff(em.table.Estimated_Marginal_Mean([4 3])'),3)),...
            ')     p = ' num2str(UFEF.pVal)])

        disp(' ')

        disp('[[ Standard vs Oddball - Within Block]]')
        disp(['UN vs EF:    ' num2str(em.table.Estimated_Marginal_Mean([2 3])'),...
            '     (' num2str(round(diff(em.table.Estimated_Marginal_Mean([2 3])'),3)),...
            ')     p = ' num2str(UNEF.pVal)])
        disp(['UF vs EN:    ' num2str(em.table.Estimated_Marginal_Mean([4 1])'),...
            '     (' num2str(round(diff(em.table.Estimated_Marginal_Mean([4 1])'),3)),...
            ')     p = ' num2str(UFEN.pVal)])

        disp('---------------------------------------------------------------------')

    end
end
end