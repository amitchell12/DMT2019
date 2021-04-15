%% A.G. Mitchell - 10.09.19
% This script takes TVA data collected for Dunhill Medical Trust funded
% study investigating visual attention and visual motor mechanisms in AD
% (McIntosh, Rossit, Mitchell, Pal & Hornberger)
% Extracts data from the TVA paradigm (previously formatted by convert_tva.r) and models each
% individual participant's data to calculated VSTM capacity (C) and processing
% speed (k)

%% Extracting data
% on mac
%dataDir = '/Users/Alex/Documents/DMT/data/formatted_TVAecc';
%anaDir = '/Users/Alex/Documents/DMT/analysis/TVA/fits';
% on desktop
dataDir = 'S:\groups\DMT\data\TVA\formatted_TVAecc';
anaDir = 'S:\groups\DMT\analysis\TVA\ecc';
cd(dataDir)
files = dir(fullfile(dataDir, 'subject*.dat'));

%% TVA file prep
for i = 1:length(files)
    cd(dataDir)
    filecheck = tvacheckdatafile(files(i).name) %leave open    
    % load data
    % STD masked conditions
    tvadata{i} = tvaloader(files(i).name, 'STD', [6 7]);
end

%% Running modelling
for i = 1:length(files)
    participantID = files(i).name(1:10);
    %% Fitting model
    [theta{i}, tvamodel{i}, tvadata{i}, df{i}] = tvafit(tvadata{i}, 8);
    %see the fitted values
    % leave tva report on data
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
end

%% Saving
cd(anaDir)
for i=1:length(files)
    if(i==1) tvalpr('OutputEcc.txt','',tvadata{i},tvamodel{i},theta{i}); 
    end
    tvalpr('OutputEcc.txt',files(i).name,tvadata{i},tvamodel{i},theta{i});
end

%% Exporting