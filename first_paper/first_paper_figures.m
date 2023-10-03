%% This script creates all figures used in my first paper


% By Andrew John Buggee



%% FIGURE 1 - Plot the ensemble MEAN of droplet size, liquid water content and 
% number concentration for non-precipitating clouds. Add an adiabatic
% profile for the liquid water content and effective radius to show the
% mean profiles are close to adiabatic, supporting my assumption


clear variables

% First, load ensemble data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    % --- non-precip profiles only, LWC>0.03, Nc>1  ----
    load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1',...
        '/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_19-Sep-2023'])


    % --- non-precip profiles only, LWC>0.03, Nc>1, stop at max LWC ----
%     load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1',...
%         '/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_stopAtMaxLWC_Nc-threshold_1_19-Sep-2023'])

    % --- non-precip profiles only, LWC>0.005, Nc>1  ----
%     load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/',...
%         'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.005_Nc-threshold_1_17-Sep-2023.mat'])

elseif strcmp(whatComputer, 'andrewbuggee')==true

    % --- non-precip profiles only, LWC>0.005, Nc>1 ----
    % load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
    %  '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.005_Nc-threshold_1_17-Sep-2023'])

    % --- non-precip profiles only, LWC>0.03, Nc>1, stop at max LWC ----
%     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
%         '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_stopAtMaxLWC_Nc-threshold_1_19-Sep-2023'])

    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
%     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
%         '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_19-Sep-2023'])

    % --- all profiles, LWC>0.005, Nc>1 ----
    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
        '/ensemble_profiles_from_14_files_LWC-threshold_0.005_Nc-threshold_1_14-Sep-2023'])


end




% using the mean ensemble function to plot the mean vertical profile of
% the ensemble



    
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Segment re, LWC, Nc into N bins along optical depth

% In order to compute a mean vertical profile, we have to first normalize
% the vertical extent so that all profiles lie between values [0,1]. Then
% we break up the vertical component in n discrete bins. Within each bin we
% can compute the mean, median and standard deviation

n_bins = 30; % number of segments the noramlized vertical component is broken up into

bin_edges = 0:1/n_bins:1;

% set up an empty cell array for all the values of each variable of interest
% within each segment boundaries. Let's do this for droplet size, total
% number concentration and liquid water content
vertically_segmented_attributes = cell(n_bins, 3);


normalized_altitude = cell(1, length(ensemble_profiles.lwc));


for nn = 1:length(ensemble_profiles.lwc)

    % first we need to normalize the vertical component of all profiles
    normalized_altitude{nn} = (ensemble_profiles.altitude{nn} - min(ensemble_profiles.altitude{nn}))./...
                               (max(ensemble_profiles.altitude{nn}) - min(ensemble_profiles.altitude{nn}));

    % the data is stored in altitude space. 

    re = ensemble_profiles.re{nn};
    lwc = ensemble_profiles.lwc{nn};
    Nc = ensemble_profiles.Nc{nn};



    % for each profile, we need to segment the variables of interest into n
    % bins.

    for bb = 1:length(bin_edges)-1

        % grab all re, LWC, and Nc values within each bin. Segment them
        % accordingly
        if bb==1
            index_segment = normalized_altitude{nn}>=bin_edges(bb) & normalized_altitude{nn}<=bin_edges(bb+1);

        else
            index_segment = normalized_altitude{nn}>bin_edges(bb) & normalized_altitude{nn}<=bin_edges(bb+1);
        end

        % store the effective radius values
        vertically_segmented_attributes{bb, 1} = [vertically_segmented_attributes{bb, 1}; re(index_segment)];

        % store the liquid water content values
        vertically_segmented_attributes{bb, 2} = [vertically_segmented_attributes{bb, 2}; lwc(index_segment)];

        % store the total droplet number concentration values
        vertically_segmented_attributes{bb, 3} = [vertically_segmented_attributes{bb, 3}; Nc(index_segment)];



    end



end



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% Create a PDF object at each level in the cloud and fit a distribution to this PDF

% store the refection of each null hypothesis and the p-value for each
% chi-squared test

re_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
re_p_normal = zeros(1, size(vertically_segmented_attributes,1));

re_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
re_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

re_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
re_p_gamma = zeros(1, size(vertically_segmented_attributes,1));


lwc_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_normal = zeros(1, size(vertically_segmented_attributes,1));

lwc_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

lwc_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_gamma = zeros(1, size(vertically_segmented_attributes,1));



Nc_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_normal = zeros(1, size(vertically_segmented_attributes,1));

Nc_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

Nc_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_gamma = zeros(1, size(vertically_segmented_attributes,1));

for bb = 1:size(vertically_segmented_attributes, 1)


    % -----------------------------------------------
    % ------- EFFECTIVE DROPLET RADIUS FITTING ------
    % -----------------------------------------------


    % fit the effective radius data to a normal distribution
    re_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'normal');
    [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb));

    % fit the effective radius data to a log-normal distribution
    re_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'lognormal');
    [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb));

    % fit the effective radius data to a gamma distribution
    re_fit_gamma(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'gamma');
    [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb));

    % make plot of normal log-normal and gamma fits
%     figure; subplot(1,3,1); plot(re_fit_normal(bb)); title('Normal Fit'); xlabel('r_e (\mum)')
%     subplot(1,3,2); plot(re_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('r_e (\mum)')
%     subplot(1,3,3); plot(re_fit_gamma(bb)); title('Gamma Fit'); xlabel('r_e (\mum)')
%     set(gcf, 'Position', [0 0 1200 500])




    % -------------------------------------------
    % ------- LIQUID WATER CONTENT FITTING ------
    % -------------------------------------------


    % fit the liquid water content data to a normal distribution
    lwc_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'normal');
    [lwc_reject_normal(bb), lwc_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF',lwc_fit_normal(bb));

    % fit the liquid water content data to a log-normal distribution
    lwc_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'lognormal');
    [lwc_reject_lognormal(bb), lwc_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF',lwc_fit_lognormal(bb));

    % fit the liquid water content data to a gamma distribution
    lwc_fit_gamma(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'gamma');
    [lwc_reject_gamma(bb), lwc_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF', lwc_fit_gamma(bb));

    % make plot of normal log-normal and gamma fits
%     figure; subplot(1,3,1); plot(lwc_fit_normal(bb)); title('Normal Fit'); xlabel('LWC (g/m^{2})')
%     subplot(1,3,2); plot(lwc_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('LWC (g/m^{2})')
%     subplot(1,3,3); plot(lwc_fit_gamma(bb)); title('Gamma Fit'); xlabel('LWC (g/m^{2})')
%     set(gcf, 'Position', [0 0 1200 500])




    % -------------------------------------------
    % ------- NUMBER CONCENTRATION FITTING ------
    % -------------------------------------------


    % fit the number concentration data to a normal distribution
    Nc_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'normal');
    [Nc_reject_normal(bb), Nc_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF',Nc_fit_normal(bb));

    % fit the number concentration content data to a log-normal distribution
    Nc_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'lognormal');
    [Nc_reject_lognormal(bb), Nc_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF',Nc_fit_lognormal(bb));

    % fit the number concentration content data to a gamma distribution
    Nc_fit_gamma(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'gamma');
    [Nc_reject_gamma(bb), Nc_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF', Nc_fit_gamma(bb));

    % make plot of normal log-normal and gamma fits
%     figure; subplot(1,3,1); plot(Nc_fit_normal(bb)); title('Normal Fit'); xlabel('N_c (cm^{-3})')
%     subplot(1,3,2); plot(Nc_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('LWC (cm^{-3})')
%     subplot(1,3,3); plot(Nc_fit_gamma(bb)); title('Gamma Fit'); xlabel('LWC (cm^{-3})')
%     set(gcf, 'Position', [0 0 1200 500])
    

end


% Now let's find the where the hypothesis was not rejected (reject_ = 0)
% which means the chi-squared test is confident in the choice of
% distribution to within 5% uncertainty

bin_names = {'Normal', 'Log-Normal', 'Gamma'};
% -----------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS FITTING ------
% -----------------------------------------------
 [max__re_p, idx_re_p] = max([re_p_normal; re_p_lognormal; re_p_gamma],[], 1);

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_re_p==1), sum(idx_re_p==2), sum(idx_re_p==3)]); 
% title('r_e best distribution fit'); ylabel('Counts')



% -------------------------------------------
% ------- LIQUID WATER CONTENT FITTING ------
% -------------------------------------------
[max__lwc_p, idx_lwc_p] = max([lwc_p_normal; lwc_p_lognormal; lwc_p_gamma],[], 1);

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_lwc_p==1), sum(idx_lwc_p==2), sum(idx_lwc_p==3)]); 
% title('LWC best distribution fit'); ylabel('Counts')


% -------------------------------------------
% ------- NUMBER CONCENTRATION FITTING ------
% -------------------------------------------

[max__Nc_p, idx_Nc_p] = max([Nc_p_normal; Nc_p_lognormal; Nc_p_gamma],[], 1);

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_Nc_p==1), sum(idx_Nc_p==2), sum(idx_Nc_p==3)]); 
% title('N_c best distribution fit'); ylabel('Counts')



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Compute the mean LWC, re, and Nc of each layer

% ---- most common best fit distribution for r_e was is the log-normal dist ---
re_logNormal_std = zeros(n_bins, 1);
re_logNormal_mean = zeros(n_bins, 1);


% ---- most common best fit distribution for LWC was is the normal dist ---
lwc_mean = zeros(n_bins, 1);
lwc_std = zeros(n_bins, 1);


% ---- most common best fit distribution for N_c was is the normal dist ---
Nc_mean = zeros(n_bins, 1);
Nc_std = zeros(n_bins, 1);

bin_center = zeros(n_bins, 1);




for bb = 1:n_bins
    

    % ----- COMPUTE STATISTICS FOR DROPLET SIZE -----

    % find the mean of the log normal distribution
    re_logNormal_mean(bb) = exp(re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2 /2);

    % find squareroot of the variance of the lognormal distribution
    re_logNormal_std(bb) = sqrt(exp(2*re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2)*(exp(re_fit_lognormal(bb).sigma^2) - 1));



    % ----- COMPUTE STATISTICS FOR LIQUID WATER CONTENT -----

    % compute the mean value for the current bin
    % the mean of the distribution (the standard way of computing the expected value)
    % is also the mean of the normal distribution. They are identical.
    lwc_mean(bb) = mean(vertically_segmented_attributes{bb,2});       % g/cm^3 - mean liqiud water content
    %lwc_mean(bb) = lwc_fit_gamma(bb).mean;       % g/cm^3 - mean liqiud water content

    % compute the standard deviation of the current bin
    % the std of the distribution (the standard way of computing the squareroot of the variance)
    % is also the std of the normal distribution. They are identical.
    lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % g/cm^3 - standard deviation
    %lwc_std(bb) = lwc_fit_gamma(bb).std;         % g/cm^3 - standard deviation



    % ----- COMPUTE STATISTICS FOR DROPLET NUMBER CONCENTRATION -----

    % compute the mean value for the current bin
    %Nc_mean(bb) = Nc_fit_gamma(bb).mean;       % cm^(-3) - mean number concentration
    Nc_mean(bb) = Nc_fit_normal(bb).mean;       % cm^(-3) - mean number concentration

    % compute the standard deviation of the current bin
    %Nc_std(bb) = Nc_fit_gamma(bb).std;         % cm^(-3) - standard deviation
    Nc_std(bb) = Nc_fit_normal(bb).std;         % cm^(-3) - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end



% ----------------------------------------------------------------------
% -------------- Make a subplot of all 3 mean profiles -----------------
% ----------------------------------------------------------------------




figure; 

% plot the mean effective radius
subplot(1,3,1)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [re_logNormal_mean - re_logNormal_std; flipud(re_logNormal_mean + re_logNormal_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(re_logNormal_mean, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter', 'latex')

% set x axis boundaries
xlim([4, 10])                   % microns



% plot the mean liquid water content profile
subplot(1,3,2)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [lwc_mean-lwc_std; flipud(lwc_mean + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(lwc_mean, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter','latex')

% set x axis boundaries
xlim([0, 0.6])                   % g/cm^3

% set the figure title
title(['Mean Profiles:  $LWC \geq$', num2str(ensemble_profiles.inputs.LWC_threshold), ' $g/m^{3}$',...
    '   $N_c \geq$',  num2str(ensemble_profiles.inputs.Nc_threshold), ' $cm^{-3}$'],...
    'Interpreter','latex')





% plot the mean droplet number concentration
subplot(1,3,3)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
x = [Nc_mean-Nc_std; flipud(Nc_mean + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(Nc_mean, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter', 'latex')

% set x axis boundaries
xlim([0, 320])                   % cm^(-3)

% set the size of the figure
set(gcf, 'Position', [0 0 1255 625])




% ----------------------------------------------------------------------
% ----------------------- Adiabatic Curve Fits -------------------------
% ----------------------------------------------------------------------

% ------------------------- LIQUID WATER CONTENT ---------------------------
% ----- Fit an adiabatic curve to the mean liquid water content profile -----
nudge_from_top =2;
lwc_slope = (lwc_mean(end-nudge_from_top) - lwc_mean(1))/(bin_center(end-nudge_from_top) - bin_center(1));
lwc_intercept = lwc_mean(1) - lwc_slope*bin_center(1);
lwc_adiabatic_fit = lwc_slope*bin_center + lwc_intercept;

% add to subplot(1,3,1)
subplot(1,3,2); hold on
plot(lwc_adiabatic_fit(1:end-nudge_from_top), bin_center(1:end-nudge_from_top), 'k', 'LineWidth',2)

% -- Include a legend in the lower right-hand corner of the 2nd subplot --
legend({'Standard Deviation', 'Mean Profile', 'Adiabatic Fit'}, 'Interpreter','latex',...
    'Location','southeast', 'FontSize', 17)



% -------------------- EFFECTIVE DROPLET RADIUS -----------------------
% Plot an adiabatic curve fit to the mean droplet radius profile 
nudge_from_bottom = 3;

% ----- Fit an adiabatic curve to the mean droplet profile -----
% use droplet profile function to create adiabatic fit
re_adiabatic_fit = create_droplet_profile2([re_logNormal_mean(end), re_logNormal_mean(1 + nudge_from_bottom)],...
                bin_center(1+nudge_from_bottom:end),'altitude', 'adiabatic');

% derive droplet profile from adiabatic fit of LWC


% add to subplot(1,3,1)
subplot(1,3,1); hold on
plot(re_adiabatic_fit, bin_center(1+nudge_from_bottom:end), 'k', 'LineWidth',2)





% ---------------------- DROPLET NUMBER CONCENTRATION ---------------------------
% ----- Fit an adiabatic curve to Number Concentration based on Adiabatic fits above -----
% Nc_adiabatic_fit = 3*(lwc_adiabatic_fit./1e6)./...
%                 (4*pi*(re_adiabatic_fit*1e-4).^3);




% --- Create a Textbox stating these are only non-precipitating clouds ---
% Create textbox
annotation('textbox',[0.134 0.802 0.142 0.114],...
    'String',{'Non-Precipitating clouds only ($LWP_{2DC}<1 \,g/m^{2}$)'},...
    'LineWidth',2,...
    'Interpreter','latex',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');





%% FIGURE 2 - Plot the ensemble MEDIAN of droplet size, liquid water content and 
% number concentration for non-precipitating clouds. Add an adiabatic
% profile for the liquid water content and effective radius to show the
% mean profiles are close to adiabatic, supporting my assumption


clear variables

% First, load ensemble data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    % --- non-precip profiles only, LWC>0.03, Nc>1  ----
    load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1',...
        '/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_19-Sep-2023'])


    % --- non-precip profiles only, LWC>0.03, Nc>1, stop at max LWC ----
%     load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1',...
%         '/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_stopAtMaxLWC_Nc-threshold_1_19-Sep-2023'])

    % --- non-precip profiles only, LWC>0.005, Nc>1  ----
%     load(['/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/',...
%         'ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.005_Nc-threshold_1_17-Sep-2023.mat'])

elseif strcmp(whatComputer, 'andrewbuggee')==true

    % --- non-precip profiles only, LWC>0.005, Nc>1 ----
    load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
     '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.005_Nc-threshold_1_17-Sep-2023'])

    % --- non-precip profiles only, LWC>0.03, Nc>1, stop at max LWC ----
%     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
%         '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_stopAtMaxLWC_Nc-threshold_1_19-Sep-2023'])

    % --- non-precip profiles only, LWC>0.03, Nc>1 ----
%     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
%         '/SPS_1/ensemble_profiles_without_precip_from_14_files_LWC-threshold_0.03_Nc-threshold_1_19-Sep-2023'])

    % --- all profiles, LWC>0.005, Nc>1 ----
%     load(['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/VOCALS_REx/vocals_rex_data',...
%         '/ensemble_profiles_from_14_files_LWC-threshold_0.005_Nc-threshold_1_14-Sep-2023'])


end




% using the median ensemble function to plot the median vertical profile of
% the ensemble



    
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Segment re, LWC, Nc into N bins along optical depth

% In order to compute a median vertical profile, we have to first normalize
% the vertical extent so that all profiles lie between values [0,1]. Then
% we break up the vertical component in n discrete bins. Within each bin we
% can compute the mean, median and standard deviation

n_bins = 30; % number of segments the noramlized vertical component is broken up into

bin_edges = 0:1/n_bins:1;

% set up an empty cell array for all the values of each variable of interest
% within each segment boundaries. Let's do this for droplet size, total
% number concentration and liquid water content
vertically_segmented_attributes = cell(n_bins, 3);


normalized_altitude = cell(1, length(ensemble_profiles.lwc));


for nn = 1:length(ensemble_profiles.lwc)

    % first we need to normalize the vertical component of all profiles
    normalized_altitude{nn} = (ensemble_profiles.altitude{nn} - min(ensemble_profiles.altitude{nn}))./...
                               (max(ensemble_profiles.altitude{nn}) - min(ensemble_profiles.altitude{nn}));

    % the data is stored in altitude space. 

    re = ensemble_profiles.re{nn};
    lwc = ensemble_profiles.lwc{nn};
    Nc = ensemble_profiles.Nc{nn};



    % for each profile, we need to segment the variables of interest into n
    % bins.

    for bb = 1:length(bin_edges)-1

        % grab all re, LWC, and Nc values within each bin. Segment them
        % accordingly
        if bb==1
            index_segment = normalized_altitude{nn}>=bin_edges(bb) & normalized_altitude{nn}<=bin_edges(bb+1);

        else
            index_segment = normalized_altitude{nn}>bin_edges(bb) & normalized_altitude{nn}<=bin_edges(bb+1);
        end

        % store the effective radius values
        vertically_segmented_attributes{bb, 1} = [vertically_segmented_attributes{bb, 1}; re(index_segment)];

        % store the liquid water content values
        vertically_segmented_attributes{bb, 2} = [vertically_segmented_attributes{bb, 2}; lwc(index_segment)];

        % store the total droplet number concentration values
        vertically_segmented_attributes{bb, 3} = [vertically_segmented_attributes{bb, 3}; Nc(index_segment)];



    end



end



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

% Create a PDF object at each level in the cloud and fit a distribution to this PDF

% store the refection of each null hypothesis and the p-value for each
% chi-squared test

re_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
re_p_normal = zeros(1, size(vertically_segmented_attributes,1));

re_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
re_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

re_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
re_p_gamma = zeros(1, size(vertically_segmented_attributes,1));


lwc_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_normal = zeros(1, size(vertically_segmented_attributes,1));

lwc_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

lwc_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
lwc_p_gamma = zeros(1, size(vertically_segmented_attributes,1));



Nc_reject_normal = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_normal = zeros(1, size(vertically_segmented_attributes,1));

Nc_reject_lognormal = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_lognormal = zeros(1, size(vertically_segmented_attributes,1));

Nc_reject_gamma = zeros(1, size(vertically_segmented_attributes,1));
Nc_p_gamma = zeros(1, size(vertically_segmented_attributes,1));

for bb = 1:size(vertically_segmented_attributes, 1)


    % -----------------------------------------------
    % ------- EFFECTIVE DROPLET RADIUS FITTING ------
    % -----------------------------------------------


    % fit the effective radius data to a normal distribution
    re_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'normal');
    [re_reject_normal(bb), re_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_normal(bb));

    % fit the effective radius data to a log-normal distribution
    re_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'lognormal');
    [re_reject_lognormal(bb), re_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF',re_fit_lognormal(bb));

    % fit the effective radius data to a gamma distribution
    re_fit_gamma(bb) = fitdist(vertically_segmented_attributes{bb,1}, 'gamma');
    [re_reject_gamma(bb), re_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,1}, 'CDF', re_fit_gamma(bb));

    % make plot of normal log-normal and gamma fits
%     figure; subplot(1,3,1); plot(re_fit_normal(bb)); title('Normal Fit'); xlabel('r_e (\mum)')
%     subplot(1,3,2); plot(re_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('r_e (\mum)')
%     subplot(1,3,3); plot(re_fit_gamma(bb)); title('Gamma Fit'); xlabel('r_e (\mum)')
%     set(gcf, 'Position', [0 0 1200 500])




    % -------------------------------------------
    % ------- LIQUID WATER CONTENT FITTING ------
    % -------------------------------------------


    % fit the liquid water content data to a normal distribution
    lwc_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'normal');
    [lwc_reject_normal(bb), lwc_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF',lwc_fit_normal(bb));

    % fit the liquid water content data to a log-normal distribution
    lwc_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'lognormal');
    [lwc_reject_lognormal(bb), lwc_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF',lwc_fit_lognormal(bb));

    % fit the liquid water content data to a gamma distribution
    lwc_fit_gamma(bb) = fitdist(vertically_segmented_attributes{bb,2}, 'gamma');
    [lwc_reject_gamma(bb), lwc_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,2}, 'CDF', lwc_fit_gamma(bb));

    % make plot of normal log-normal and gamma fits
%     figure; subplot(1,3,1); plot(lwc_fit_normal(bb)); title('Normal Fit'); xlabel('LWC (g/m^{2})')
%     subplot(1,3,2); plot(lwc_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('LWC (g/m^{2})')
%     subplot(1,3,3); plot(lwc_fit_gamma(bb)); title('Gamma Fit'); xlabel('LWC (g/m^{2})')
%     set(gcf, 'Position', [0 0 1200 500])




    % -------------------------------------------
    % ------- NUMBER CONCENTRATION FITTING ------
    % -------------------------------------------


    % fit the number concentration data to a normal distribution
    Nc_fit_normal(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'normal');
    [Nc_reject_normal(bb), Nc_p_normal(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF',Nc_fit_normal(bb));

    % fit the number concentration content data to a log-normal distribution
    Nc_fit_lognormal(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'lognormal');
    [Nc_reject_lognormal(bb), Nc_p_lognormal(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF',Nc_fit_lognormal(bb));

    % fit the number concentration content data to a gamma distribution
    Nc_fit_gamma(bb) = fitdist(vertically_segmented_attributes{bb,3}, 'gamma');
    [Nc_reject_gamma(bb), Nc_p_gamma(bb)] = chi2gof(vertically_segmented_attributes{bb,3}, 'CDF', Nc_fit_gamma(bb));

    % make plot of normal log-normal and gamma fits
%     figure; subplot(1,3,1); plot(Nc_fit_normal(bb)); title('Normal Fit'); xlabel('N_c (cm^{-3})')
%     subplot(1,3,2); plot(Nc_fit_lognormal(bb)); title('Log-Normal Fit'); xlabel('LWC (cm^{-3})')
%     subplot(1,3,3); plot(Nc_fit_gamma(bb)); title('Gamma Fit'); xlabel('LWC (cm^{-3})')
%     set(gcf, 'Position', [0 0 1200 500])
    

end


% Now let's find the where the hypothesis was not rejected (reject_ = 0)
% which means the chi-squared test is confident in the choice of
% distribution to within 5% uncertainty

bin_names = {'Normal', 'Log-Normal', 'Gamma'};
% -----------------------------------------------
% ------- EFFECTIVE DROPLET RADIUS FITTING ------
% -----------------------------------------------
 [max__re_p, idx_re_p] = max([re_p_normal; re_p_lognormal; re_p_gamma],[], 1);

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_re_p==1), sum(idx_re_p==2), sum(idx_re_p==3)]); 
% title('r_e best distribution fit'); ylabel('Counts')



% -------------------------------------------
% ------- LIQUID WATER CONTENT FITTING ------
% -------------------------------------------
[max__lwc_p, idx_lwc_p] = max([lwc_p_normal; lwc_p_lognormal; lwc_p_gamma],[], 1);

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_lwc_p==1), sum(idx_lwc_p==2), sum(idx_lwc_p==3)]); 
% title('LWC best distribution fit'); ylabel('Counts')


% -------------------------------------------
% ------- NUMBER CONCENTRATION FITTING ------
% -------------------------------------------

[max__Nc_p, idx_Nc_p] = max([Nc_p_normal; Nc_p_lognormal; Nc_p_gamma],[], 1);

% figure; histogram('Categories', bin_names, 'BinCounts', [sum(idx_Nc_p==1), sum(idx_Nc_p==2), sum(idx_Nc_p==3)]); 
% title('N_c best distribution fit'); ylabel('Counts')



% ----------------------------------------------------------------------
% ----------------------------------------------------------------------
% ----------------------------------------------------------------------

%  Compute the Median LWC, re, and Nc of each layer

% ---- most common best fit distribution for r_e was is the log-normal dist ---
re_logNormal_std = zeros(n_bins, 1);
re_logNormal_median = zeros(n_bins, 1);
re_custom_logNormal_std_larger = zeros(n_bins, 1);
re_custom_logNormal_std_smaller = zeros(n_bins, 1);

% ---- most common best fit distribution for LWC was is the normal dist ---
lwc_median = zeros(n_bins, 1);
lwc_std = zeros(n_bins, 1);


% ---- most common best fit distribution for N_c was is the normal dist ---
Nc_median = zeros(n_bins, 1);
Nc_std = zeros(n_bins, 1);

bin_center = zeros(n_bins, 1);




for bb = 1:n_bins
    

    % ----- COMPUTE STATISTICS FOR DROPLET SIZE -----

    % find the mean of the log normal distribution
    re_logNormal_median(bb) = re_fit_lognormal(bb).median;

    % find squareroot of the variance of the lognormal distribution
    re_logNormal_std(bb) = sqrt(exp(2*re_fit_lognormal(bb).mu + re_fit_lognormal(bb).sigma^2)*(exp(re_fit_lognormal(bb).sigma^2) - 1));

    % Let's also compute the average deviation from the median value when
    % radii are larger and smaller than the median.
    % For radii larger than the median...
    idx_larger = vertically_segmented_attributes{bb,1}>re_logNormal_median(bb);
    re_larger = vertically_segmented_attributes{bb,1}(idx_larger);
    re_custom_logNormal_std_larger(bb) = sqrt(mean((re_larger - re_logNormal_median(bb)).^2));
    %re_custom_logNormal_std_larger(bb) = mean(re_larger - re_logNormal_median(bb));

    % For radii smaller than the median...
    idx_smaller = vertically_segmented_attributes{bb,1}<re_logNormal_median(bb);
    re_smaller = vertically_segmented_attributes{bb,1}(idx_smaller);
    re_custom_logNormal_std_smaller(bb) = sqrt(mean((re_smaller - re_logNormal_median(bb)).^2));
    %re_custom_logNormal_std_smaller(bb) = mean(re_logNormal_median(bb) - re_smaller);



    % ----- COMPUTE STATISTICS FOR LIQUID WATER CONTENT -----

    % compute the mean value for the current bin
    % the mean of the distribution (the standard way of computing the expected value)
    % is also the mean of the normal distribution. They are identical.
    lwc_median(bb) = median(vertically_segmented_attributes{bb,2});       % g/cm^3 - mean liqiud water content

    % compute the standard deviation of the current bin
    % the std of the distribution (the standard way of computing the squareroot of the variance)
    % is also the std of the normal distribution. They are identical.
    lwc_std(bb) = std(vertically_segmented_attributes{bb,2});         % g/cm^3 - standard deviation
    %lwc_std(bb) = lwc_fit_gamma(bb).std;         % g/cm^3 - standard deviation



    % ----- COMPUTE STATISTICS FOR DROPLET NUMBER CONCENTRATION -----

    % compute the mean value for the current bin
    Nc_median(bb) = Nc_fit_normal(bb).median;       % cm^(-3) - mean number concentration

    % compute the standard deviation of the current bin
    %Nc_std(bb) = Nc_fit_gamma(bb).std;         % cm^(-3) - standard deviation
    Nc_std(bb) = Nc_fit_normal(bb).std;         % cm^(-3) - standard deviation

    % compute the bin center which is the tau location of the mean data
    bin_center(bb) = (bin_edges(bb+1) - bin_edges(bb))/2 + bin_edges(bb);

end



% ----------------------------------------------------------------------
% -------------- Make a subplot of all 3 mean profiles -----------------
% ----------------------------------------------------------------------




figure; 

% plot the median effective radius
subplot(1,3,1)

% plot the standard deviation of the mean profile as an transparent area
% centered around the mean radius profile
%x = [re_logNormal_median - re_logNormal_std; flipud(re_logNormal_median + re_logNormal_std)];
x = [re_logNormal_median - re_custom_logNormal_std_smaller; flipud(re_logNormal_median + re_custom_logNormal_std_larger)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the mean droplet profile
plot(re_logNormal_median, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<r_e(z)>$ $(\mu m)$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter', 'latex')

% set x axis boundaries
xlim([4, 12])                   % microns



% plot the median liquid water content profile
subplot(1,3,2)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
x = [lwc_median-lwc_std; flipud(lwc_median + lwc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(lwc_median, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<LWC(z)>$ $(g/m^{3})$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter','latex')

% set x axis boundaries
xlim([0, 0.6])                   % g/cm^3

% set the figure title
title(['Median Profiles:  $LWC \geq$', num2str(ensemble_profiles.inputs.LWC_threshold), ' $g/m^{3}$',...
    '   $N_c \geq$',  num2str(ensemble_profiles.inputs.Nc_threshold), ' $cm^{-3}$'],...
    'Interpreter','latex')





% plot the median droplet number concentration
subplot(1,3,3)

% plot the standard deviation of the median profile as an transparent area
% centered around the mean radius profile
x = [Nc_median-Nc_std; flipud(Nc_median + Nc_std)];
y = [bin_center; flipud(bin_center)];
fill(x,y,mySavedColors(2,'fixed'), 'EdgeAlpha', 0, 'FaceAlpha', 0.2)

hold on

% plot the median droplet profile
plot(Nc_median, bin_center, 'Color', mySavedColors(2, 'fixed'))


grid on; grid minor
xlabel('$<N_c(z)>$ $(cm^{-3})$', 'Interpreter','latex')
ylabel('Normalized Altitude', 'Interpreter', 'latex')

% set x axis boundaries
xlim([0, 320])                   % cm^(-3)

% set the size of the figure
set(gcf, 'Position', [0 0 1255 625])




% ----------------------------------------------------------------------
% ----------------------- Adiabatic Curve Fits -------------------------
% ----------------------------------------------------------------------

% ------------------------- LIQUID WATER CONTENT ---------------------------
% ----- Fit an adiabatic curve to the mean liquid water content profile -----
nudge_from_top = 3;
lwc_slope = (lwc_median(end-nudge_from_top) - lwc_median(1))/(bin_center(end-nudge_from_top) - bin_center(1));
lwc_intercept = lwc_median(1) - lwc_slope*bin_center(1);
lwc_adiabatic_fit = lwc_slope*bin_center + lwc_intercept;

% add to subplot(1,3,1)
subplot(1,3,2); hold on
plot(lwc_adiabatic_fit(1:end-nudge_from_top), bin_center(1:end-nudge_from_top), 'k', 'LineWidth',2)

% -- Include a legend in the lower right-hand corner of the 2nd subplot --
legend({'Standard Deviation', 'Median Profile', 'Adiabatic Fit'}, 'Interpreter','latex',...
    'Location','southeast', 'FontSize', 17)



% -------------------- EFFECTIVE DROPLET RADIUS -----------------------
% Plot an adiabatic curve fit to the mean droplet radius profile 
nudge_from_bottom = 2;

% ----- Fit an adiabatic curve to the mean droplet profile -----
% use droplet profile function to create adiabatic fit
re_adiabatic_fit = create_droplet_profile2([re_logNormal_median(end), re_logNormal_median(1 + nudge_from_bottom)],...
                bin_center(1+nudge_from_bottom:end),'altitude', 'adiabatic');

% derive droplet profile from adiabatic fit of LWC


% add to subplot(1,3,1)
subplot(1,3,1); hold on
plot(re_adiabatic_fit, bin_center(1+nudge_from_bottom:end), 'k', 'LineWidth',2)





% ---------------------- DROPLET NUMBER CONCENTRATION ---------------------------
% ----- Fit an adiabatic curve to Number Concentration based on Adiabatic fits above -----
% density = 1;            % g/cm^3
% Nc_adiabatic_fit = 3*(lwc_adiabatic_fit./1e6)./...
%                 (4*pi*density*(re_adiabatic_fit*1e-4).^3);          % #/cm^3



subplot(1,3,3)
hold on;
%plot(Nc_adiabatic_fit, bin_center, 'k', 'LineWidth', 2)
xline(mean(Nc_median), 'Color', 'black', 'LineWidth', 2, 'Alpha',1,...
    'Label', [num2str(round(mean(Nc_median))), ' $cm^{-3}$'], 'Interpreter','latex',...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','top',...
    'Fontsize', 20)



% --- Create a Textbox stating these are only non-precipitating clouds ---
% Create textbox
% annotation('textbox',[0.134 0.802 0.07 0.06],...
%     'String',{'All Profiles'},...
%     'LineWidth',2,...
%     'Interpreter','latex',...
%     'FontSize',17,...
%     'FontName','Helvetica Neue',...
%     'FitBoxToText','off');


annotation('textbox',[0.134 0.802 0.142 0.114],...
    'String',{'Non-Precipitating clouds only ($LWP_{2DC}<1 \,g/m^{2}$)'},...
    'LineWidth',2,...
    'Interpreter','latex',...
    'FontSize',17,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');




%% FIGURE 3

% Plot weighting functions of the first 7 MODIS spectral channels
% These weighting functions were created using the VOCALS-REx data set from
% Nov-9-2023. The droplet profile and optical depth were modeled after the
% vertical profile sampled at 1.734 hrs after the plane took off. The SZA
% was set as the value measured by MODIS for the median pixel, the pixel
% found closest the C130 aircraft in the middle of it's ascent through the
% cloud. 

filenames = {'2D_MC_05-Sep-2023_Wavelength_469_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
             '2D_MC_05-Sep-2023_Wavelength_555_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
             '2D_MC_05-Sep-2023_Wavelength_645_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
             '2D_MC_05-Sep-2023_Wavelength_858.5_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
             '2D_MC_05-Sep-2023_Wavelength_1240_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
             '2D_MC_05-Sep-2023_Wavelength_1640_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat',...
             '2D_MC_05-Sep-2023_Wavelength_2130_N-Photons_10000000_N-Layers_100_Tau0_15_SZA_27.mat'};

% plot as a pdf
probability_str = 'pdf';

% define the wavelengths as the changing variables
wavelength = modisBands(1:7);
changing_variable = wavelength(:,1);        % center wavelenghts in nm

plot_probability_scat_top_maxDepth_with_changing_variable(filenames, probability_str ,changing_variable)



%% FIGURE 3

% Plot the single wavelength retrieved droplet radius for each of the 7 wavelengths

figure; 
plot(inputs.dropletProfile.re, inputs.dropletProfile.tau_layer_mid_points)
% flip the y axis so 0 is at the depth
set(gca, 'YDir', 'reverse')
grid on; grid minor
xlabel('$r_e(\tau)$','Interpreter','latex');
ylabel('$\tau$','Interpreter','latex')

hold on
for nn = 1:length(filenames)
    
    % Load a simulation
    load(filenames{nn})

    % compute the PDF of photons scattering out the top after reaching a
    % max depth of tau
        % First select those photons that were scattered out the top

    index_scatter_out_top = final_state.scatter_out_top_INDEX;

    [scatter_out_top_maxDepth_PDF, scatter_out_top_maxDepth_PDF_edges] = ...
        histcounts(photon_tracking.maxDepth(index_scatter_out_top),'Normalization',probability_str);

    % independent TAU variable for the PDF
    tau_pdf = scatter_out_top_maxDepth_PDF_edges(1:end-1) + diff(scatter_out_top_maxDepth_PDF_edges)/2;

    % interpolate to get a droplet profile the same lenght as the PDF
    interp_re = interp1(inputs.dropletProfile.tau_layer_mid_points, inputs.dropletProfile.re, tau_pdf, 'spline');

    retrieved_re = trapz(tau_pdf, interp_re .* scatter_out_top_maxDepth_PDF);

    % add a vertical line to the plot
    xline(retrieved_re, 'LineStyle','--', 'LineWidth',2, 'Color', 'k',...
        'Label',['$\lambda = $',  num2str(round(inputs.mie.wavelength(1))), ' $\mu m$'],...
        'LabelHorizontalAlignment','center', 'LabelVerticalAlignment','bottom', ...
        'Interpreter','latex', 'FontSize', 15)



end

set(gcf, 'Position',[0 0 1000 630])


%% FIGURE 4 

% Compute the standard deviation of droplet size over some identified
% length scale. Slide this length scale across the horizontal profile and
% compute the standard deviation for each window. Do this for each profile
% in a data set


clear variables

% Plot the liquid wtaer content, effective radius and number concentration
% as a function of horizontal distance travelled for a single day of data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    folder_path = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/';

elseif strcmp(whatComputer, 'andrewbuggee')==true

    folder_path = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/', ...
    'VOCALS_REx/vocals_rex_data/SPS_1/'];


end

% Oct-15-2008 Data
%filename = 'RF01.20081015.164800_201200.PNI.nc';

% Oct-18-2008 Data
%filename = 'RF02.20081018.130300_213000.PNI.nc';

% Oct-21-2008 Data
%filename = 'RF03.20081021.060000_142400.PNI.nc';

% Oct-25-2008 Data
% filename = 'RF05.20081025.062900_152500.PNI.nc';

% ----- November 9 data -----
filename = 'RF11.20081109.125700_213600.PNI.nc';

% ------ November 11 data -----
%filename = 'RF12.20081111.125000_214500.PNI.nc';

% ----- November 13 data -----
%filename = 'RF13.20081113.125700_215700.PNI.nc';

% ----- November 15 data -----
%filename = 'RF14.20081115.125800_220400.PNI.nc';


vocalsRex = readVocalsRex([folder_path, filename]);


% ---- set thresholds for the LWC and Nc ---
LWC_threshold = 0.03;       % g/m^3
Nc_threshold = 1;           % cm^{-3}
max_vertical_displacement = 10;     % meters

% ---- Find all Horizontal Profiles ---
horz_profs = find_horizontalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, Nc_threshold, max_vertical_displacement);

% --- Define the horizontal length scale to compute statistics over ---
length_scale = 1000;        % meters

% Loop through each profile

for nn = 1:length(horz_profs.lwc)

    % step trough each individual droplet size profile
    mean_lengthScale{nn} = [];
    std_lengthScale{nn} = [];

    for rr = 2:length(horz_profs.re{nn})

        % First, check to see if the remaining distance between rr and the
        % end of our horizontal profile is greater than our length scale. 
        % If not, break the for loop
        if (horz_profs.horz_dist{nn}(end) - horz_profs.horz_dist{nn}(rr))>length_scale
        
            idx_displacement = 0;
            % find data points that make up 1 km
            dist_1km = horz_profs.horz_dist{nn}(rr+idx_displacement) - horz_profs.horz_dist{nn}(rr-1);
    
            while dist_1km<length_scale
                
                % step to the next data point
                idx_displacement = idx_displacement + 1;
                dist_1km = horz_profs.horz_dist{nn}(rr+idx_displacement) - horz_profs.horz_dist{nn}(rr-1);
    
            end
    
            % when the distance between two data points reaches 1000 meters,
            % stop and calculate the mean and standard deviation
            mean_lengthScale{nn} = [mean_lengthScale{nn}, mean(horz_profs.re{nn}(rr-1 : rr+idx_displacement))];         % microns
            std_lengthScale{nn} = [std_lengthScale{nn}, std(horz_profs.re{nn}(rr-1 : rr+idx_displacement))];           % microns

        else

            break

        end

    end

end



% ----- PLOT 1 -----
% Plot the mean versus the standard deviation for each profile
figure; 
for nn = 1:length(mean_lengthScale)
    plot(std_lengthScale{nn}, mean_lengthScale{nn}, '.-', 'Linewidth', 1.75,...
        'MarkerSize', 15)
    hold on

end
    
    
    % Include an x axis label on the middle plot
xlabel('$\sigma(r_e)$ ($\mu m$)', 'Interpreter','latex');

grid on; grid minor;
ylabel('$\left<r_e\right>$ ($\mu m$)', 'Interpreter','latex')

title(['Mean and STD over length scale of ', num2str(length_scale), ' meters'], 'Interpreter','latex')

% set plot size
set(gcf, 'Position', [0 0 1200 625])



% ----- PLOT 2 -----
% Try plotting the standard deviation for every length scale as a histogram
legend_str = cell(1, length(std_lengthScale));

figure;
for nn = 1:length(mean_lengthScale)
    
    histogram(std_lengthScale{nn}, 10, 'FaceAlpha', 0.5)
    hold on

    legend_str{nn} = ['index = ', num2str(nn), ' $\left<\sigma(r_e)\right> = \,$', num2str(mean(std_lengthScale{nn})), ' $\mu m$'];



end

% Include an x axis label on the middle plot
xlabel('$\sigma(r_e)$ ($\mu m$)', 'Interpreter','latex');

grid on; grid minor;
ylabel('Counts', 'Interpreter','latex')

title(['STD over length scale of ', num2str(length_scale), ' meters'], 'Interpreter','latex')

legend(legend_str, 'Interpreter','latex', 'Location','best', 'FontSize', 19)

% set plot size
set(gcf, 'Position', [0 0 1200 625])




%% FIGURE 5

clear variables

% Plot the liquid wtaer content, effective radius and number concentration
% as a function of horizontal distance travelled for a single day of data

% grab the filepath name according to which computer is being used

if strcmp(whatComputer, 'anbu8374')==true

    folder_path = '/Users/anbu8374/Documents/MATLAB/HyperSpectral_Cloud_Retrieval/VOCALS_REx/vocals_rex_data/SPS_1/';

elseif strcmp(whatComputer, 'andrewbuggee')==true

    folder_path = ['/Users/andrewbuggee/Documents/MATLAB/CU Boulder/Hyperspectral_Cloud_Retrievals/', ...
    'VOCALS_REx/vocals_rex_data/SPS_1/'];


end

% Oct-15-2008 Data
%filename = 'RF01.20081015.164800_201200.PNI.nc';

% Oct-18-2008 Data
%filename = 'RF02.20081018.130300_213000.PNI.nc';

% Oct-21-2008 Data
%filename = 'RF03.20081021.060000_142400.PNI.nc';

% Oct-25-2008 Data
% filename = 'RF05.20081025.062900_152500.PNI.nc';

% ----- November 9 data -----
filename = 'RF11.20081109.125700_213600.PNI.nc';

% ------ November 11 data -----
%filename = 'RF12.20081111.125000_214500.PNI.nc';

% ----- November 13 data -----
%filename = 'RF13.20081113.125700_215700.PNI.nc';

% ----- November 15 data -----
%filename = 'RF14.20081115.125800_220400.PNI.nc';


vocalsRex = readVocalsRex([folder_path, filename]);


% ---- set thresholds for the LWC and Nc ---
LWC_threshold = 0.03;       % g/m^3
Nc_threshold = 1;           % cm^{-3}
max_vertical_displacement = 10;     % meters

% ---- Find all Horizontal Profiles ---
horz_profs = find_horizontalProfiles_VOCALS_REx(vocalsRex, LWC_threshold, Nc_threshold, max_vertical_displacement);


% ---- Plot the properties of all horizontal profiles collected ---
plot_horiztonal_profiles_LWC_and_re_and_Nc(horz_profs, 1:length(horz_profs.re), false)


