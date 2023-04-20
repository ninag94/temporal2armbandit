%% main protocol


% choose the directory containing the files to analyse
Directory = uigetdir;
prefix = input('provide the prefix of the files to select:'); %f.e. 2_NosePoke
AllSessionFilesList = dir(fullfile(Directory,[prefix,'*mat']));
prompt2 = "Enter the path of the directory for saving the plots:";   % f.e. '/Users/ninagrimme/Desktop/Matlab/#2_OneArmCued/plots'
path = input(prompt2);

%maybe ask user for a path for storing the results f.e. plots etc.
%(if-statement in case the user does not want to store the results)

%define mode for either analyzing session per session or all sessions
mode = input('Provide the type of analysis: single session or all sessions: [s|a]');
if mode ~= 's' && mode ~= 'a'
    mode = 's';
    errormsg = 'Invalid type of analysis. The default setting is a single session analysis.';
    disp(errormsg)
end

% Initialize an empty cell array to store the data
AllSessionFiles = cell(1, numel(AllSessionFilesList));

% Loop through each file and load the data into the cell array
for i = 1:numel(AllSessionFilesList)
    file_name = fullfile(Directory, AllSessionFilesList(i).name);
    %check the size of the file
    if AllSessionFilesList(i).bytes < 100000
        continue
    end
   
    AllSessionFiles{i} = load(file_name);
end

% Convert the cell array to a numeric array
indxemptycells = any(~cellfun('isempty',AllSessionFiles), 1);
AllSessionFiles = AllSessionFiles(indxemptycells);
AllSessionFiles = cell2mat(AllSessionFiles);


%switch statement for mode = s (session per session analysis)
switch mode
    case 's'
        for i = 1:length(AllSessionFiles)
            File = AllSessionFiles(i);
            DateSession = File.SessionData.Custom.General.SessionDate;
            Animal = File.SessionData.Custom.General.Subject;
            % all functions for analyzing the data go here!
            %[A,B] = notbaitedanalysis(File); % A, B is currently only a placeholder

            %running average plot
            [MChoiceLeftSmoothed,RunningAveragePlot]= blockanalysis(File);
            Plot = fullfile(path,sprintf('%s_%s_running_average.png',Animal,DateSession));
            saveas(RunningAveragePlot,Plot);

            %all session events plot
            [AllSessionEvents,AllEventsPlot] = alltrials_fixedwithholding(File);
            Plot = fullfile(path,sprintf('%s_%s_allsessionevents.png',Animal,DateSession));
            saveas(AllEventsPlot,Plot);
           

            
        end
  

    % and mode = a (all sessions pooled)
    case 'a'
        %get the date of first and last session for the period of time for
        %naming plots etc.
        DateSessionFirst = AllSessionFiles(1).SessionData.Custom.General.SessionDate;
        DateSessionLast = AllSessionFiles(numel(AllSessionFiles)).SessionData.Custom.General.SessionDate;
        Animal = AllSessionFiles(1).SessionData.Custom.General.Subject;

        %all functions for analyzing the data go here

        %all functions for plotting go here


end





%call temporal2armbandit_dataanalysis.mat




%call temporal2armbandit_plot.mat



%clean up working directory