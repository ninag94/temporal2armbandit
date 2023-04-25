%% 2ArmBandittask analysis function

function FigHandle = Analysis(DataFile)
global BpodSystem

if nargin < 1
    DataFile = BpodSystem.Data;
end

%% Load related data to local variabels
try
    Animal = str2double(DataFile.Info.Subject);
catch
    Animal = -1;
end

nTrials = DataFile.nTrials;
ChoiceLeft = DataFile.Custom.TrialData.ChoiceLeft;
Baited = DataFile.Custom.TrialData.Baited;
IncorrectChoice = DataFile.Custom.TrialData.IncorrectChoice;
NoChoice = DataFile.Custom.TrialData.NoDecision;
NoTrialStart = DataFile.Custom.TrialData.NoTrialStart;
BrokeFix = DataFile.Custom.TrialData.BrokeFixation;
EarlyWith = DataFile.Custom.TrialData.EarlyWithdrawal;
SkippedFeedback = DataFile.Custom.TrialData.SkippedFeedback;
Rewarded = DataFile.Custom.TrialData.Rewarded;

%FeedbackWaitingTime = DataFile.Custom.TrialData.FeedbackWaitingTime;
FeedbackWaitingTime = rand(684,1)*10'; %delete this
FeedbackWaitingTime = FeedbackWaitingTime';  %delete this
RewardProb = DataFile.Custom.TrialData.RewardProb;
LightLeft = DataFile.Custom.TrialData.LightLeft;

ChoiceLeftRight = [ChoiceLeft; 1-ChoiceLeft]; 
%% Plot based on task design/ risk type
switch DataFile.SettingsFile.GUIMeta.RiskType.String{DataFile.SettingsFile.GUI.RiskType}
    case 'Fix'
    %%
        %not yet implemented

    case 'BlockRand'
    %%
        %not yet implemented

    case 'BlockFix'
    %%
        %not yet implemented

    case 'BlockFixHolding'
    %%
        FigHandle = figure('Position', [360  187	1056	598],...
                           'NumberTitle', 'off',...
                           'Name', num2str(Animal));

        %running choice average

        RewardProbLeft = RewardProb(1,:);
        
        %% Block switching behaviour across trial
        subplot(2,2,1) %needs adjustment! depends on the number of plots in the end
        if ~isempty(ChoiceLeft) && ~all(isnan(ChoiceLeft))
            xdata = 1:nTrials;
            plot(xdata, RewardProbLeft, '-k', 'Color', [.5,.5,.5], 'LineWidth', 2);
            hold on;

            ChoiceLeftSmoothed = smooth(ChoiceLeft, 10, 'moving','omitnan'); %current bin width: 10 trials
            plot(xdata, ChoiceLeftSmoothed, '-k', 'LineWidth', 2);
            ylim([0 1]);
            xlim([0 nTrials]);
            ylabel('Ratio of Left Choices (%)')
            xlabel('Trials')
            title('Block switching behviour')
        end

        %% all trials overview (counts for each observed behavior)

        ChoiceLeftRight = [ChoiceLeft; 1-ChoiceLeft];
        NotBaited = any((Baited == 0) .* ChoiceLeftRight, 1);

        counts = zeros(1,6)';
        events = {'NoChoice', 'BrokeFix','EarlyWith','SkippedFeedback','Rewarded','NotBaited'};
        AllSessionEvents = table(counts,'RowNames',events);

        if ~isempty(NoChoice) && ~all(isnan(NoChoice))
            AllSessionEvents('NoChoice','counts') = {length(NoChoice(NoChoice==1))};
        end
        if ~isempty(BrokeFix) && ~all(isnan(BrokeFix))
            AllSessionEvents('BrokeFix','counts') = {length(BrokeFix(BrokeFix==1))};
        end
        if ~isempty(EarlyWith) && ~all(isnan(EarlyWith))
            AllSessionEvents('EarlyWith','counts') = {length(EarlyWith(EarlyWith==1))};
        end
        if ~isempty(SkippedFeedback) && all(isnan(SkippedFeedback))
            AllSessionEvents('SkippedFeedback','counts') = {length(SkippedFeedback(SkippedFeedback==1))}; %-length(indxNotBaited(indxNotBaited==1))};
        end                                                                                              %is skipped feedback now seperate from not-baited?
        if ~isempty(Rewarded) && ~all(isnan(Rewarded))
            AllSessionEvents('Rewarded','counts') = {length(Rewarded(Rewarded==1))};
        end
        if ~isempty(NotBaited) && ~all(isnan(NotBaited))

            AllSessionEvents('NotBaited','counts')= {length(NotBaited(NotBaited==1))};

        end
    
    
        y = AllSessionEvents.counts;  
        xlabels = {'No Choice','Broke Fix','Early With','Skipped Feedback','Rewarded','NotBaited'};
        colors = {'yellow', "#77AC30",'blue',"#EDB120",'black','cyan'};  

        subplot(2,2,2)   %needs adjustment, depends on the number of subplots
        hold on

        for i = 1:length(y)
            bar(i, y(i), 'FaceColor', colors{i}, 'EdgeColor', 'none')
        end

        set(gca, 'XTick', 1:length(y))
        set(gca, 'XTickLabel', xlabels)
        %legend(xlabels)
        title('All trials');
        ylabel("counts");



        %% distribution of waiting time of notbaited trials
        
        if ~isempty(FeedbackWaitingTime) && ~all(isnan(FeedbackWaitingTime))
            
            RewardProbChosen = RewardProb .* ChoiceLeftRight;
            RewardProbChosen = RewardProbChosen(1,:) + RewardProbChosen(2,:);
            WaitingTimeBlocks = [FeedbackWaitingTime;RewardProbChosen];
            NotBaitedWT = WaitingTimeBlocks(:,NotBaited==1);
            PHigh = max(RewardProbChosen);
            PLow = min(RewardProbChosen);
            WTHigh = NotBaitedWT(:,NotBaitedWT(2,:) == PHigh);
            WTLow = NotBaitedWT(:,NotBaitedWT(2,:) == PLow);
            nHigh = length(WTHigh);
            nLow = length(WTLow);
            meanHigh = [mean(WTHigh(1));PHigh];
            meanLow = [mean(WTLow(1));PLow];

            subplot(2,2,3);    %needs adjustment!

            scatter(WTLow(2,:),WTLow(1,:),'cyan');
            hold on
            scatter(WTHigh(2,:),WTHigh(1,:),'blue');
            plot(meanHigh(2,:),meanHigh(1,:),'x','MarkerSize',10,'MarkerEdgeColor','black','LineWidth',2);
            plot(meanLow(2,:),meanHigh(2,:),'x','MarkerSize',10,'MarkerEdgeColor','black','LineWidth',2);
            xlim([0 1]);
            ticks = [PLow,PHigh];
            xticks(ticks);
            xticklabels({'PLow','PHigh'});
            ylabel('time investment (s)');
            text1 = sprintf('n = %d',nLow);
            text2 = sprintf('n = %d',nHigh);
            text(PLow,max(FeedbackWaitingTime),text1);
            text(PHigh,max(FeedbackWaitingTime),text2);
            title('Waiting times for not-baited trials')
            

        end

        %% calculating the drinking times

        DrinkingTime = [];

        for i = 1:nTrials

            if Rewarded(i) == 1

                DrinkingTime(end+1) = DataFile.RawEvents.Trial{i}.States.Drinking(2) - DataFile.RawEvents.Trial{1,1}.States.Drinking(1);

            end

        end

        if ~all(isnan(DrinkingTime))

            subplot(2,2,4);              % specification of the binning could be added
            histogram(DrinkingTime,'FaceColor',[.5,.5,.5],'EdgeColor',[1,1,1]);
            xlabel('drinking times (s)')
            ylabel('n')
            title('Distribution of drinking times')

        end


    case 'Cued'

        FigHandle = figure('Position',[ 360         187        1056         598],'NumberTitle','off','Name',Animal);

        %waiting time of not-baited trials per reward probability

        %get only the actually used reward probabilities

        LightLeftRight = [LightLeft;1-LightLeft];
        LightRewardProb = RewardProb .* LightLeftRight;
        RewardProbUsed = LightRewardProb(1,:) + LightRewardProb(2,:);

        if ~isempty(FeedbackWaitingTime)

            %get the waiting time and reward probability for not-baited trials

            indxNotBaited = (IncorrectChoice==0) & any((Baited == 0) .* ChoiceLeftRight, 1); 
            NotBaitedWT = FeedbackWaitingTime(indxNotBaited==1); 
            NotBaitedRP = RewardProbUsed(indxNotBaited==1);
            WaitingTime = NotBaitedWT';
            RewardProb = NotBaitedRP';
            NotBaitedWTRP = table(RewardProb,WaitingTime); % does not need to be in this format, can also be an array/matrix
            WTRPSummary = groupsummary(NotBaitedWTRP,"RewardProb",{"mean","std"}); %then groupstats need to be used instead of groupsummary

            %plot waiting time and reward probability

            subplot(2,2,1)                      %needs adjustment to the final number of plots
            xdata = WTRPSummary.RewardProb;
            ydata = WTRPSummary.mean_WaitingTime;
            swarmchart(NotBaitedWTRP.RewardProb,NotBaitedWTRP.WaitingTime,'cyan','filled')
            hold on
            errbar = WTRPSummary.std_WaitingTime;
            plot(xdata,ydata,'diamond','MarkerSize',8,'MarkerEdgeColor','black','LineWidth',2);
            title("Not Baited trials")
            xlabel("reward probability in %");
            ylabel("time investment in s ");
            errorbar(xdata,ydata,errbar,'Linestyle','none','Color','k');
            xdataticks = WTRPSummary.RewardProb';
            xticks(xdataticks)           
            xticklabels(string(xdataticks));        %test this?
            xlim([0 1])
            counts = string(WTRPSummary.GroupCount);
            ylabels = ydata + 0.2;   
            text(xdata,ylabels,counts,'Color','black','FontSize',12)
            annotation('textbox', [0.77, 0.7, 0.1, 0.1], 'String', "total counts",'Color','black')

        end
end