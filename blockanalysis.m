%% function to analyse the data of the decision behavior within and between the different blocks
%how does the change of reward prob affect behavior
%also i-1 trial effect maybe in this function?

%% function to analyse the choices in relation to the blocks and switching between blocks
% matching behavior?

%for the later plot: x-axis: trial number
%plot 1: reawrd probability as a line (makes blocks very clear)
%plot 2: choices for which reward probability did the animal decide in each
%trial?

%function [y,yx] = blockanalysis(x) not yet implemented

%gathering the information that is neeeded
Trialnumber = 1:SessionData.nTrials;
RewardProb = SessionData.Custom.TrialData.RewardProb; %should be sufficient to plot the blocks
ChoiceLeft = SessionData.Custom.TrialData.ChoiceLeft;

RewardProbLeft = RewardProb(1,:);
ChoiceLeftRight = [ChoiceLeft; 1-ChoiceLeft]; 
ChoiceRewardProb = RewardProb .* ChoiceLeftRight;
%ChoiceRewardProb = [ChoiceRewardProb(1,:) + ChoiceRewardProb(2,:);Trialnumber];
%ChoiceRewardProb = ChoiceRewardProb(:,~isnan(ChoiceRewardProb(1,:))); %keep only the trials with a decision


p = figure;

plot(Trialnumber,RewardProbLeft,'-k','Color',[.5,.5,.5],'LineWidth',2);
ylim([0 1]);
hold on
%plot(ChoiceRewardProb(2,:),ChoiceRewardProb(1,:));


