%% A.G. Mitchell - 10.09.19
% This script takes TVA data collected for Dunhill Medical Trust funded
% study investigating visual attention and visual motor mechanisms in AD
% (McIntosh, Rossit, Mitchell, Pal & Hornberger)
% Extracts data from the TVA paradigm (previously formatted by convert_tva.r) and models each
% individual participant's data to calculated VSTM capacity (C) and processing
% speed (k)

%% Extracting data
% on mac
%dataDir = '/Users/Alex/Documents/DMT/data/formatted_TVA';
%anaDir = '/Users/Alex/Documents/DMT/analysis/TVA/fits';
% on desktop
dataDir = 'S:\groups\DMT\data\TVA\formatted_TVA';
anaDir = 'S:\groups\DMT\analysis\TVA\all';
cd(dataDir)
files = dir(fullfile(dataDir, 'subject*.dat'));

%% TVA file prep
for i = 1:length(files)
   
    filecheck = tvacheckdatafile(files(i).name) %leave open
    
    % load data
    % need to include the 'STD' [5 6] for masked conditions
    tvadata{i} = tvaloader(files(i).name, 'STD', [6 7]);
end
    
    %% Fitting model
for i = 1:length(files)
     participantID = files(i).name(1:10);
    [theta{i}, tvamodel{i}, tvadata{i}, df{i}] = tvafit(tvadata{i}, [], 'FREE');
    %see the fitted values
    tvareport(tvadata{i}, tvamodel{i}, theta{i});
    
    %% Plotting
    [tt, oo, pp, cc] = tvaplot(tvadata{i}, tvamodel{i}, theta{i});
    % first figure - by condition (unmasked at end)
    figure()
    plot(cc,pp,'*')
    ax = gca;
    xLabels = tt;
    xticklabels(ax, xLabels); 
    ylabel('Mean score (k)'); xlabel('Exposure duration (ms)');
    xlim([0 8]); ylim([0 4]);
    
    cd([anaDir filesep 'fit-plots'])
    png_name = sprintf('%s_%s.png', participantID);
    saveas(gcf, png_name)
    
    % second figure by exposure duration (unmasked intermixed)
%     figure()
%     plot(tt,pp,'*')
%     ylabel('Mean score (k)'); xlabel('Exposure duration (ms)');
%     xlim([0 250]); ylim([0 4]);
    
    %% Saving individual data
    %cd([anaDir filesep 'fits'])
    
    %datafilename = sprintf('%s_%s_fits.txt', participantID, filename);
    %tvaexport(datafilename, 'DIR', 'fit.mat');
end

cd(anaDir)
for i=1:length(files)
    if(i==1) tvalpr('Output.txt','',tvadata{i},tvamodel{i},theta{i}); 
    end
    tvalpr('Output.txt',filename(i).name,tvadata{i},tvamodel{i},theta{i});
end

%% Exporting