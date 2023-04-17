%% function to analyse the data/create the data output for the plotting of the notbaited trials

function [NotBaitedWTRP,StatsWTRP] = notbaitedanalysis(x) %probably multiple outputs for the diff plots
    %get all the trial information
    Baited = x.SessionData.Custom.TrialData.Baited;
    IncorrectChoice = x.SessionData.Custom.TrialData.IncorrectChoice;
    FeedbackWaitingTime = x.SessionData.Custom.TrialData.FeedbackWaitingTime;
    RewardProb = x.SessionData.Custom.TrialData.RewardProb;
    LightLeft = x.SessionData.Custom.TrialData.LightLeft;
    ChoiceLeft = x.SessionData.Custom.TrialData.ChoiceLeft;

    % get only the actually used reward probabilities
    LightLeftRight = [LightLeft;1-LightLeft];
    LightRewardProb = RewardProb .* LightLeftRight;
    RewardProbUsed = LightRewardProb(1,:) + LightRewardProb(2,:);

    %% relationship between not-baited trials and reward probability?
    % get the indices of not-baited trials
    ChoiceLeftRight = [ChoiceLeft; 1-ChoiceLeft]; 
    indxNotBaited = (IncorrectChoice~=1) & any((Baited == 0) .* ChoiceLeftRight, 1); 
    % get the waiting time for non-baited trials
    NotBaitedWT = FeedbackWaitingTime(indxNotBaited==1); 
    % get the reward probability for the non-baited trials
    NotBaitedRP = RewardProbUsed(indxNotBaited==1);
    NotBaitedWTRP = [NotBaitedWT;NotBaitedRP]; % here I'm not so sure whether to use a table or a matrix
    [mean,std,count,grps] = grpstats(NotBaitedWTRP(:,1), NotBaitedWTRP(:,2), {'mean', 'std', 'numel','gname'});
    StatsWTRP = [grps;mean;std;count]; %is this actually in the right order?

