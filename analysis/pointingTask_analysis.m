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
%nParticipants = 2:4; %number of participants
nParticipants = 2; %for testing
% task 1 = closed loop beep, task 2 = closed loop fix, task 3 = closed loop
% norm, task 4 = open loop
nTasks = 1:4; %number of tasks used during testing

allData = struct; %data structure

%% Extracting data
cd(dataPath) %to where the data is
for p = 1:length(nParticipants)
    ppID = sprintf('hp%0*d',2,nParticipants(p)); %for use when navigating files
    anaFileName = sprintf('%s_pointingAnalysis.mat',ppID); %name of final analysis file
    
    dirData = [dataPath filesep ppID]; %where to extract data
    % making analysis folder for data extraction
    dirAna = mkdir([anaPath filesep ppID]); 
   
    % importing relevant files
    cd(dirData)
    fileName = sprintf('subject90%d_pointingTask_*.csv', nParticipants(p));
    d = dir(fileName);
    names = {d.name}; % directory with important filenames
    
    % Loop through each task importing and analysing critical data
    for t = 1:length(nTasks)
        %task name - generally done alphabetically
        switch t
            case 1
                taskName = 'CLbeep'; %take name
            case 2
                taskName = 'CLfix';
            case 3
                taskName = 'CLnorm';
            case 4
                taskName = 'OL';
        end
        %importing data
        [~,~,left_data] = xlsread(names{t*2-1});
        [~,~,right_data] = xlsread(names{t*2});  
        % sorting so everything is in the right order
        [temp, order] = sort(left_data(1,:)); %left data
        left_data = left_data(:,order);
        [temp, order] = sort(right_data(1,:)); %right data
        right_data = right_data(:,order);

        % adding to structure
        allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).left = left_data;
        allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).right = right_data;
        
        % Extracting important variables from the data
        %eyeMove = 

        if t == 2
            %land_x = 
        else
        end
       
    end
end