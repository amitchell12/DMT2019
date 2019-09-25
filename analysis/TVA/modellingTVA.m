%% A.G. Mitchell - 10.09.19
% This script takes TVA data collected for Dunhill Medical Trust funded
% study investigating visual attention and visual motor mechanisms in AD
% (McIntosh, Rossit, Mitchell, Pal & Hornberger)
% Extracts data from the TVA paradigm (previously formatted by convert_tva.r) and models each
% individual participant's data to calculated VSTM capacity (C) and processing
% speed (k)

%% Extracting data
dataDir = 'M:\Alex_Files\Experiments\DMT2019\raw_data\formatted_TVA';
anaDir = 'M:\Alex_Files\Experiments\DMT2019\analysed_data\';
cd(dataDir)
files = dir(fullfile(dataDir, 'subject*.dat'));

%% TVA file prep
for i = 1:length(files)
    cd(dataDir)
    participantID = files(i).name(1:10);
    filecheck = tvacheckdatafile(files(i).name) %leave open
    % load data
    % need to include the 'STD' [5 6] for masked conditions
    tvadata = tvaloader(files(i).name, 'STD', [5 6]);
    % tva report on data
    tvareport(tvadata) %leave open to view
    
    %% Fitting model
    [theta, tvamodel, tvadata, df] = tvafit(tvadata);
    %see the fitted values
    tvareport(tvadata, tvamodel, theta);
    
    %% Plotting
    [tt, oo, pp, cc] = tvaplot(tvadata, tvamodel, theta);
    % first figure - by condition (unmasked at end)
    figure()
    plot(cc,pp,'*')
    ax = gca;
    xLabels = tt;
    xticklabels(ax, xLabels); 
    ylabel('Mean score (k)'); xlabel('Exposure duration (ms)');
    xlim([0 8]); ylim([0 4]);
    
    % second figure by exposure duration (unmasked intermixed)
    figure()
    plot(tt,pp,'*')
    ylabel('Mean score (k)'); xlabel('Exposure duration (ms)');
    xlim([0 250]); ylim([0 4]);
    
    %% Saving individual data
    cd(anaDir)
    save fit theta tvamodel tvadata
    datafilename = sprintf('%s_fits.csv', participantID);
    tvaexport(datafilename, 'DIR', 'fit.mat');
end

%% Exporting