%all trials 2armbandit

function [AllSessionEvents, p] = alltrials_fixedwithholding(x)

 nTrials = x.SessionData.nTrials;
 ChoiceLeft = x.SessionData.Custom.TrialData.ChoiceLeft;
 Baited = x.SessionData.Custom.TrialData.Baited;
 NoChoice = x.SessionData.Custom.TrialData.NoDecision;
 NoTrialStart = x.SessionData.Custom.TrialData.NoTrialStart;
 BrokeFix = x.SessionData.Custom.TrialData.BrokeFixation;
 EarlyWith = x.SessionData.Custom.TrialData.EarlyWithdrawal;
 SkippedFeedback = x.SessionData.Custom.TrialData.SkippedFeedback;
 Rewarded = x.SessionData.Custom.TrialData.Rewarded;


 ChoiceLeftRight = [ChoiceLeft; 1-ChoiceLeft];
 indxNotBaited = any((Baited == 0) .* ChoiceLeftRight, 1);

 counts = zeros(1,6)';
 events = {'NoChoice', 'BrokeFix','EarlyWith','SkippedFeedback','Rewarded','NotBaited'};
 AllSessionEvents = table(counts,'RowNames',events);

 if ~isempty(NoChoice)
     AllSessionEvents('NoChoice','counts')= {length(NoChoice(NoChoice==1))};
 end
 if ~isempty(BrokeFix)
     AllSessionEvents('BrokeFix','counts')= {length(BrokeFix(BrokeFix==1))};
 end
 if ~isempty(EarlyWith)
     AllSessionEvents('EarlyWith','counts')= {length(EarlyWith(EarlyWith==1))};
 end
 if ~isempty(SkippedFeedback)
     AllSessionEvents('SkippedFeedback','counts')= {length(SkippedFeedback(SkippedFeedback==1))}; %-length(indxNotBaited(indxNotBaited==1))};
 end                                                                                              %is skipped feedback now seperate from not-baited?
 if ~isempty(Rewarded)  
     AllSessionEvents('Rewarded','counts')= {length(Rewarded(Rewarded==1))};
 end
 if ~isempty(indxNotBaited)
     AllSessionEvents('NotBaited','counts')= {length(indxNotBaited(indxNotBaited==1))};
 end
    
    
y = AllSessionEvents.counts;  
xlabels = {'No Choice','Broke Fix','Early With','Skipped Feedback','Rewarded','NotBaited'};
colors = {'yellow', "#77AC30",'blue',"#EDB120",'black','cyan'};   
p = figure;
hold on

for i = 1:length(y)
    bar(i, y(i), 'FaceColor', colors{i}, 'EdgeColor', 'none')
end

set(gca, 'XTick', 1:length(y))
set(gca, 'XTickLabel', xlabels)
%legend(xlabels)
title([Animal,"-",DateSession],"FontSize",12);
ylabel("counts");