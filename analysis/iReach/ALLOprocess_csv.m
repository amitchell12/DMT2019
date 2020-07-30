% allocentric iReach ipad task processing code for csv files  
% 15/7/20
% check trajL variable to make sure real movements being selected not just
% noise (will have v small values I assume if noise)

close all; clear;
trialPlot=0; %% will do trial by trial plots of speed X time, will be much slower. Set to 0 to not do it and 1 to do it.
trajPlot=0;%% will do trial by trial plots of x by y, will be much slower. Set to 0 to not do it and 1 to do it.
sampleThreshold=5; %% some subjects e.g. 213 have only one sample for a given trial - cut trials with < 5 
%(will appear in lower left at zero,zero in end point accuracy plot)

% get fileNames from the user here by getting them to click on touch and
% other file
[FileName,PathName,FilterIndex] = uigetfile('.csv','***Pick _ALLO_TOUCH file***'); % user load main filename
[FileName2,PathName2,FilterIndex2] = uigetfile('.csv', '***Pick the other ALLO file***'); % user load second filename

% read ipad csv files
[NUM,TXT,RAW]=xlsread([PathName FileName]); % breaks it into numeric, text and raw data
[NUM2,TXT2,RAW2]=xlsread([PathName2 FileName2]);
fp=find(NUM2(:,5)==1); %% last trial labelled 1 taken as first main trial
fp=fp(end);
reactionTime=NUM2(fp:end,6); %% already explicit for main trials 
subject_code=num2str(NUM(1,1));
targets=[103.89091 65.11818; 88.89091 25.11818; 43.89091 65.11818; 58.89091 25.11818]; %target positions

%%%%%%%% change coordinates here for different subjects %%%%%%%
endpoints=[38.89091 55.11818; 53.89091 95.11818; 93.89091 95.11818; 108.89091 55.11818]; %perfect end-points for above

% first step - ignore practice trials!
a=strfind(TXT(:,15),'false');% find real trials
for i=1:length(a)
    
    s(i)=isempty(a{i}); %% find onset of real trial
    
    if(s(i)==0)
        break;
    end
        
end

data=NUM(i-1:end,:); %puts real data into another matrix

nTrials=max(data(:,5)); % trial numbers
cutInd=[];

% loop through each trial
for trial=1:nTrials
    
    % extract data per trial
    f=find(data(:,5)==trial); %finds all the data for one trial
    
    if(length(f)<sampleThreshold)
        cutInd=[cutInd trial];
        continue; %% cut this trial - poor number of samples
    end
    
    trialData=data(f,:);% all the data for one trial saved into trialData
    
    hh=isnan(trialData(:,6)); %% cut the NaN's (these are from the written label in time column....)
    trialData=trialData(~hh,:);
    
    if(strcmp(subject_code,'201') || strcmp(subject_code,'202') || strcmp(subject_code,'203'))
        % solves the issue that on some trials the time index starts
        % extrememly high (3000) before going to zero
        hh=find(trialData(:,6)==0);
        trialData=trialData(hh:end,:);
    end
    
    % compute velocity criteria for start and end at each time step
    dXdY=diff(trialData(:,8:9)); %% difference in x / y position
    mag=sqrt(sum(dXdY.^2,2));  %% magnitude
    dt=diff(trialData(:,6))./1000; %% time differences in secs (presume msec originally!!!!)
    
    s=mag./dt; %% velocity (mm per sec)
    trialData(2:end,17)= s>50;  %(50mm per sec threshold for onset / offset)
    trialData(2:end,18)=s; %% store actual velocity
    
    f2=find(trialData(:,17)==1); %% movements > threshold
    if(isempty(f2))
        cutInd=[cutInd trial];
        continue; % if no movement exceeds threshold, skip to next trial
    end
    
    moveData(trial,1:2)=trialData(f2(1),8:9);  %% start position
    
    % check on end position as sometimes single movements way later in
    % file
    tmp=trialData(:,17); 
    dd=diff(tmp);
    ff=find(dd~=0);
    try
%         if(~isnan(tmp(1)))
             endInd=ff(2); % if there is some noise in the movement
%         else
%             endInd=ff(3);
%         end
    catch
        endInd=length(dd); % else take the last index
    end
    
    % store movement onset and offset, target position and trajectory for
    % later plotting etc
    moveData(trial,3:4)=trialData(endInd,8:9); %% end position
    [moveData(trial,5),IndPeak]=max(s(f2(1):endInd)); %% peak Velocity
    
    
    targetPos(trial,1:4)=trialData(f2(1),10:13); %% target position
    
    traj{trial}=[trialData(f2(1):endInd,[6 8 9 10 11 12 13])];  %% store trajectories
    trajL(trial)=length(traj{trial}); % lenght of trajectory in frames if too small then its probably an exclude (lifted pen from tablet)
    
    % time to peak velocity
    moveData(trial,6)=traj{trial}(IndPeak,1)-traj{trial}(1,1); %% time to peak velocity
    
    % calculate movement Time (calculated from travel time)
    moveT(trial)=trialData(endInd,6)-trialData(f2(1),6);
    
    % plot speed by time
    if(trialPlot)
        figure, plot(traj{trial}(:,1),s(f2(1):endInd),'b*');
        axis([0 1000 0 600])
        title(sprintf('Velocity by time_%s', num2str(trial)))

    end
    
    % plot XY trajectory
    if(trajPlot)
        figure, plot(traj{trial}(:,2),traj{trial}(:,3),'ko-'); % plot trajectories per trial
        axis([0 200 0 250])
        title(sprintf('XY Trajectory_%s', num2str(trial)))
    end
    
end  %% end the loop over trials

ff=find(cutInd==40); % if trial 40 is cut, only 39 trials will be in the output, otherwise trials marked as zeros!!!!
if(ff)
    reactionTime(40)=[];
    
end

% calculate accuracy?
targetPosCorr=[];
for i=1:4
    
    f=find(round(targetPos(:,1),2)==round(targets(i,1),2));
    targetPosCorr(f,1:2)=repmat(endpoints(i,:),length(f),1);
    
end
Acc=targetPosCorr(:,1:2)-moveData(:,3:4); %% target - end position
AE=sum(sqrt(Acc(:,1:2).^2),2); % one value per trial (absolute error)

% plot 
%xx=unique(targetPos,'rows');
figure, scatter(moveData(:,3),moveData(:,4)); %% end position
hold, plot(endpoints(:,1),endpoints(:,2),'b*') % markers for end-point positions
title('End Position by perfect end-positions')
print('-dtiff',sprintf('%s_EndPosition.tiff',subject_code));

figure, hist(moveT), title('Histogram of Movement Times');  %% histogram of movement times
print('-dtiff',sprintf('%s_MovementTime.tiff',subject_code));

figure, hist(reactionTime), title('Histogram of Reaction Times'); 
print('-dtiff',sprintf('%s_ReactionTime.tiff',subject_code));

figure, scatter(targetPos(:,1),AE), title('Absolute error by x target position');  %% scatter of AE
print('-dtiff',sprintf('%s_Absolute error.tiff',subject_code));

% save things?
outname=[subject_code '_' FileName2(1:end-4) '_processed.mat'];
outData=[targetPos targetPosCorr reactionTime moveT' trajL' moveData(:,3:6) Acc AE];
%% xtarget, Ytarget, ref x, ref y, correctendpointx, correctendpointy, RT, MT, trajL, x,y, PV, TPV,  xerror, yerror, AE

xlswrite([subject_code '_allo_processed.xlsx'],outData);

save(outname,'outData','moveData','traj','trajL','targetPos','moveT','reactionTime','nTrials','Acc','AE','cutInd');
