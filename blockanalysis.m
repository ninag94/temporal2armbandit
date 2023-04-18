
%% function to analyse the choices in relation to the blocks and switching between blocks


function [ChoiceLeftSmoothed,p] = blockanalysis(x)

%gathering the information that is neeeded

Trialnumber = 1:x.SessionData.nTrials;
RewardProb = x.SessionData.Custom.TrialData.RewardProb;
ChoiceLeft = x.SessionData.Custom.TrialData.ChoiceLeft;

RewardProbLeft = RewardProb(1,:);

%plot the data: running choice average
p = figure;

plot(Trialnumber,RewardProbLeft,'-k','Color',[.5,.5,.5],'LineWidth',2);
ylim([0 1]);
xlim([0 max(Trialnumber)]);
hold on
ChoiceLeftSmoothed = smooth(ChoiceLeft, 10, 'moving','omitnan'); %current binning: 10 trials
plot(Trialnumber,ChoiceLeftSmoothed,'-k','LineWidth',2);
ylabel('Ratio of Left Choices)')
xlabel('Trials')
hold off


