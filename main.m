%% main protocol


% set up working directory

%define mode for either analyzing session per session or all sessions

input(mode, 'Provide the type of analysis: single session or all sessions: [s|a]')
if mode goes wrong
    mode = 's';
    errormsg = 'Invalid type of analysis. The default setting is a single session analysis.';
    disp(errormsg)
end


%call temporal2armbandit_dataanalysis.mat




%call temporal2armbandit_plot.mat



%clean up working directory