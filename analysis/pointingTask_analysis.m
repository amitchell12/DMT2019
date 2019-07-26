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
% screen information
% getting pixels per cm/deg
sd = 40;
xRes = 1960; yRes = 1080;
x = 31; y = 17.5;
pix_perDeg = 43.76;
deg_perPix = 1/pix_perDeg;

% caluclating cm per pixel
pixPer_cm_x = xRes/x;
pixPer_cm_y = yRes/y;
pixPer_cm = (pixPer_cm_x+pixPer_cm_y)/2;
cm_perPix = 1/pixPer_cm;

%deg_perPix = atan2(cm_perPix/(2*sd));

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
                
        %importing data
        [~,~,left_data] = xlsread(names{t*2-1});
        [~,~,right_data] = xlsread(names{t*2});  
        % sorting so everything is in the right order
        [temp, order] = sort(left_data(1,:)); %left data
        left_data = left_data(:,order);
        [temp, order] = sort(right_data(1,:)); %right data
        right_data = right_data(:,order);
        
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

        % adding to structure
        allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).left = left_data;
        allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).right = right_data;
        
        switch t 
            case 1
                % Extracting important variables from the data
                %for the LHS
                %eyeMove = cell2mat(left_data(2:end,8)); %eye movement data - for filtering 
                target_location = cell2mat([left_data(2:end,14), left_data(2:end,15)]); %actual target location        
                touch_response = cell2mat([left_data(2:end,5), left_data(2:end,6)]); %landing position of the point to target, x and y
                trial = [1:length(target_location)]';
                leftMat = [trial, target_location, touch_response];
                allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).leftMat = leftMat;
                
                %for the RHS
                %eyeMove = cell2mat(right_data(2:end,8)); %eye movement data - for filtering 
                target_location = cell2mat([right_data(2:end,14), right_data(2:end,15)]); %actual target location        
                touch_response = cell2mat([right_data(2:end,5), right_data(2:end,6)]); %landing position of the point to target, x and y
                trial = [1:length(target_location)]';
                rightMat = [trial, target_location, touch_response];
                allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).rightMat = rightMat;
        
            case 2
                %for the LHS
                %eyeMove = cell2mat(left_data(2:end,5)); %eye movement data - for filtering 
                target_location = cell2mat([left_data(2:end,14), left_data(2:end,15)]); %actual target location        
                touch_response = cell2mat([left_data(2:end,7), left_data(2:end,8)]); %landing position of the point to target, x and y
                trial = [1:length(target_location)]';
                leftMat = [trial, target_location, touch_response];
                allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).leftMat = leftMat;
                
                %for the RHS
                %eyeMove = cell2mat(right_data(2:end,5)); %eye movement data - for filtering 
                target_location = cell2mat([right_data(2:end,14), right_data(2:end,15)]); %actual target location        
                touch_response = cell2mat([right_data(2:end,7), right_data(2:end,8)]); %landing position of the point to target, x and y
                trial = [1:length(target_location)]';
                rightMat = [trial, target_location, touch_response];
                allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).rightMat = rightMat;
                
            case 3
                %for the LHS
                %eyeMove = cell2mat(left_data(2:end,5)); %eye movement data - for filtering 
                target_location = cell2mat([left_data(2:end,13), left_data(2:end,14)]); %actual target location        
                touch_response = cell2mat([left_data(2:end,4), left_data(2:end,5)]); %landing position of the point to target, x and y
                trial = [1:length(target_location)]';
                leftMat = [trial, target_location, touch_response];
                allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).leftMat = leftMat;
                
                %for the RHS
                %eyeMove = cell2mat(right_data(2:end,5)); %eye movement data - for filtering 
                target_location = cell2mat([right_data(2:end,13), right_data(2:end,14)]); %actual target location        
                touch_response = cell2mat([right_data(2:end,4), right_data(2:end,5)]); %landing position of the point to target, x and y
                trial = [1:length(target_location)]';
                rightMat = [trial, target_location, touch_response];
                allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).rightMat = rightMat;
                
            case 4
                %for the LHS
                %eyeMove = cell2mat(left_data(2:end,5)); %eye movement data - for filtering 
                target_location = cell2mat([left_data(2:end,14), left_data(2:end,15)]); %actual target location        
                touch_response = cell2mat([left_data(2:end,5), left_data(2:end,6)]); %landing position of the point to target, x and y
                trial = [1:length(target_location)]';
                leftMat = [trial, target_location, touch_response];
                allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).leftMat = leftMat;
                
                %for the RHS
                %eyeMove = cell2mat(right_data(2:end,5)); %eye movement data - for filtering 
                target_location = cell2mat([right_data(2:end,14), right_data(2:end,15)]); %actual target location        
                touch_response = cell2mat([right_data(2:end,5), right_data(2:end,6)]); %landing position of the point to target, x and y
                trial = [1:length(target_location)]';
                rightMat = [trial, target_location, touch_response];
                allData.(sprintf('%s', ppID)).(sprintf('%s', taskName)).rightMat = rightMat;
        end
        
        %% analysing the data
        % getting error relative to target midpoint
        % CL beep, left hand data
        leftx_error = allData.(sprintf('%s', ppID)).CLbeep.leftMat(:,2) - ...
            allData.(sprintf('%s', ppID)).CLbeep.leftMat(:,4); 
        leftx_error_cm = leftx_error*cm_perPix; %converting from pix to cm
        leftx_error_deg = leftx_error*deg_perPix; %converting from pix to deg
        lefty_error = allData.(sprintf('%s', ppID)).CLbeep.leftMat(:,5) - ...
            allData.(sprintf('%s', ppID)).CLbeep.leftMat(:,3); 
        lefty_error_cm = lefty_error*cm_perPix;
        lefty_error_deg = lefty_error*deg_perPix; %converting from pix to deg
        
        % right hand data
        rightx_error = allData.(sprintf('%s', ppID)).CLbeep.rightMat(:,4) - ...
            allData.(sprintf('%s', ppID)).CLbeep.rightMat(:,2); 
        rightx_error_cm = rightx_error*cm_perPix; %converting from pix to cm
        rightx_error_deg = rightx_error*deg_perPix; %converting from pix to deg
        righty_error = allData.(sprintf('%s', ppID)).CLbeep.rightMat(:,5) - ...
            allData.(sprintf('%s', ppID)).CLbeep.rightMat(:,3); 
        righty_error_cm = righty_error*cm_perPix;
        righty_error_deg = righty_error*deg_perPix; %converting from pix to deg
        
        % saving to matrices in allData
        allData.(sprintf('%s', ppID)).CLbeep.leftError(:,1) = leftx_error_deg;
        allData.(sprintf('%s', ppID)).CLbeep.leftError(:,2) = lefty_error_deg;
        allData.(sprintf('%s', ppID)).CLbeep.rightError(:,1) = rightx_error_deg;
        allData.(sprintf('%s', ppID)).CLbeep.rightError(:,2) = righty_error_deg;
        
        
    end
end