%% AG. Mitchell - 25.07.19
%% DMT2019 Pointing task analysis script
% Takes pointing data from .csv files extracted from opensesame experiment
% Experiment: peripheral pointing task, touch targets presented at
% 24,31,37deg eccentricity
% Subtracts touch response to target from actual target location
% Does this for each target location (3 trials per target loc), for each
% eccentricity (9 trials per) and for each side (L/R - 27 trials per)

% Task will be used to compare pointing in AD patients with control
% patients
% For now- use for healthy pilot data

%% Data files and variables
dataPath = 'M:\Alex_Files\Experiments\DMT2019\DMT2019_rawData';
anaPath = 'M:\Alex_Files\Experiments\DMT2019\dataAnalysis';
nParticipants = 1:4; %number of participants
% task 1 = closed loop norm, task 2 = closed loop fix, task 3 = closed loop
% beep, task 4 = open loop
nTasks = 1:4; %number of tasks used during testing

allData = struct; %data structure

%% Extracting data
cd(dataPath) %to where the data is
for p = 1:length(nParticipants)
    ppID = sprintf('hp%0*d',2,nParticipants(p)); %for use when navigating files
    anaFileName = sprintf('%s_pointingAnalysis.mat',ppID); %name of final analysis file
    
    dirData = [dataPath filesep ppID]; %where to extract data
    % Data file names for each task
    for t = 1:length(nTasks)
        switch t
            case 1
                dataFileName_left = sprintf('subject90%d_pointingTask_CLnorm_left*', t);
            case 2
            case 3
            case 4
        end
    end
    % making analysis folder for data extraction
    dirAna = mkdir([anaPath filesep ppID]); 
   
    % importing relevant filenames
    
end