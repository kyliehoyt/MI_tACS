clear
%close all
clc
%% Load subject data
Subject_12 = load('Subj012.mat', 'sub').sub;
Subject_16 = load('Subj016.mat', 'sub').sub;
Subject_17 = load('Subj017.mat', 'sub').sub;
fs = Subject_12.Pre.restingState.run.header.SampleRate;
load("ch32Locations.mat");


% Subject groups
% tRNS_Subjects = Subject_12;
% tACS_Subjects = [Subject_16, Subject_17];
All_Subjects = [Subject_12, Subject_16, Subject_17];
All_Subjects = cleanSubjects(All_Subjects, fs, 0, [1 50]);

% Global variables
global chan_map
chan_map = string(Subject_12.Pre.restingState.run.header.Label(1:32));
alpha = [8 12];

%% Are there changes in the alpha band power before and after stimulation?
% Compute changes in alpha power before and after stimulation
norm_rs_alpha = buildPowerTable(All_Subjects, fs, {'tRNS', 'tACS', 'tACS'}, alpha);
% Determine if any of the changes are statistically significant
[ttest_results, signif_chans] = rsPowerttest(norm_rs_alpha);

%% Visualization
% plot the power change of the top changing channels for tRNS vs tACS
plot_chan = map_chan(["C3", "C4", "FZ", "FC2", "CP5", "CP6"]);
createPowerBarPlot(norm_rs_alpha, plot_chan);
% plot the power change of all channels on a topoplot for tRNS vs tACS
createTopoplots(norm_rs_alpha, ch32Locations);


%% Method using pwelch and t-test

function rs_power_change = buildPowerTable(subjects, fs, treatment_labels, band)
% Computes the total alpha band power then normalizes it to each subject's
% pre-stimulation alpha band power
% Inputs: 
%   subjects          = 1xn_sub struct of subject structs
%   fs                = sampling frequency (Hz)
%   treatment_labels  = 1xn_sub cell array of 'labels' of treatment groups in the same order as subjects
%   ch32Locations     = EEG locations .mat file
% Outputs:
%   rs_power_change   = struct with change in band power (%) for each channel and subject 
%   (values) and treatment labels (treatment)
    n_sub = length(subjects); 
    sessions = {'Pre', 'Post'};
    n_sess = length(sessions);
    n_chan = 32;
    [time, treatment] = deal(cell(n_sub*n_sess, 1)); % session labels: session time and treatment type
    session_band_powers = zeros(n_sub*n_sess, n_chan); % total alpha powers at each channel (col) for each session (row) 
    for sub = 1:n_sub
        for sess = 1:n_sess
            session = sessions{sess}; % get session time name: 'Pre' or 'Post'
            period = getfield(subjects, {sub}, session, 'restingState', 'run', 'eeg'); % resting state run
            period = trimRun(period, n_chan); % remove trailing zeros and extra channels if present
            [pxx, ~] = pwelch(period(10*fs:end-10*fs,:), fs, fs/2, 1:50, fs); % 1-50 Hz psd of middle 40 seconds at 1 Hz resolution
            total_band_power = sum(pxx(band(1):band(2), :), 1); % sum over 8-12 Hz: total alpha band power at each channel (col)
            element = (sub-1)*n_sess + sess; % next index
            session_band_powers(element, :) = total_band_power;
            treatment{element, 1} = treatment_labels{sub};
            time{element, 1} = session;
        end
    end
    rs_power = struct('values', session_band_powers, 'treatment', {treatment}, 'time', {time});
    rs_power_change = normalizePower(rs_power); % normalize to pre of each subject
end

function [ttest_results, signif_chans] = rsPowerttest(power_changes)
% Performs 2 sample t-tests for independent means to determine if there is
% a difference in the power change between tRNS subjects and tACS subjects
% at each channel
% Inputs: 
%   power_changes = struct with change in resting state alpha power (%) 
%   for each channel and subject (values), treatment labels (treatment), 
%   and session time labels (time)
% Outputs:
%   ttest_results    = p-values for each channel
%   signif_chans    = struct with names (chan), channel numbers (chan_num),
%   and p-values (p) of channels that were significantly different in total
%   power change
    n_chan = 32;
    ttest_results = zeros(n_chan, 1);
    tACSselect = cellfun(@(m)isequal(m,"tACS"),power_changes.treatment);
    tRNSselect = cellfun(@(m)isequal(m,"tRNS"),power_changes.treatment);
    for chan = 1:n_chan
        [~, ttest_results(chan)] = ttest2(power_changes.values(tRNSselect, chan), power_changes.values(tACSselect, chan));
    end
    signif_chans.chan = map_chan(ttest_results < 0.05);
    signif_chans.chan_num = find(ttest_results < 0.05);
    signif_chans.p = ttest_results(ttest_results < 0.05);
end

function bandpower_change = normalizePower(bandpowers)
% Converts the total band power of each session to the percent change
% between the pre and post stimulation sessions for each subject
% Inputs: 
%   bandpowers = struct with total band power for each channel and subject 
%   (values), treatment labels (treatment), and session time labels (time)
% Outputs:
%   bandpower_change = struct with change in band power (%) for each channel and subject 
%   (values) and treatment labels (treatment)
    Pre = cellfun(@(m)isequal(m,"Pre"),bandpowers.time);
    Post = cellfun(@(m)isequal(m,"Post"),bandpowers.time);
    pre_powers = bandpowers.values(Pre, :);
    bandpower_change.values = (bandpowers.values(Post,:)-bandpowers.values(Pre,:))./pre_powers.*100;
    bandpower_change.treatment = bandpowers.treatment(Post);
end


function createPowerBarPlot(power_changes, chan_num)
% Plots the change in band power (%) for tRNS subjects vs tACS subjects for
% the specified channels
% Inputs:
%   power_changes = struct with change in band power (%) for each channel and subject 
%   (values) and treatment labels (treatment)
%   chan_num = index of channel(s) to plot
    tACSselect = cellfun(@(m)isequal(m,"tACS"),power_changes.treatment);
    tRNSselect = cellfun(@(m)isequal(m,"tRNS"),power_changes.treatment);
    tRNS = power_changes.values(tRNSselect, chan_num);
    tACS = power_changes.values(tACSselect, chan_num);
    means = [mean(tRNS, 1); mean(tACS, 1)];
    sems = [std(tRNS, 0, 1)/sqrt(size(tRNS, 1)); std(tACS, 0, 1)/sqrt(size(tACS, 1))];
    xlab = categorical({'tRNS', 'tACS'});
    xlab = reordercats(xlab,{'tRNS', 'tACS'});
    h = bar(xlab, means);
    hold on
    for ch = 1:length(chan_num)
        errorbar(h(ch).XEndPoints, means(:, ch), sems(:, ch),'LineStyle','none','Color','k','LineWidth',2)
    end
    legend(h, map_chan(chan_num), 'FontSize',13, 'location', 'best');
    ylabel("Alpha Power Change (%)", 'FontSize', 13)
    ax = gca;
    ax.YGrid = "on";
    ax.YMinorGrid = "on";
    ax.FontSize = 11;
end

function createTopoplots(power_changes, ch32Locations)
% Plots the change in band power (%) for tRNS subjects vs tACS subjects for
% all channels on a topoplot
% Inputs:
%   power_changes = struct with change in band power (%) for each channel and subject 
%   (values) and treatment labels (treatment)
%   ch32Locations = channel locations .mat file
    tRNSselect = cellfun(@(m)isequal(m,"tRNS"),power_changes.treatment);
    tACSselect = cellfun(@(m)isequal(m,"tACS"),power_changes.treatment);
    maplim = [-20 120];
    figure();
    t = tiledlayout(1,2, 'TileSpacing','compact');
    nexttile;
    tRNS = mean(power_changes.values(tRNSselect, :), 1);
    topoplot(tRNS, ch32Locations, 'maplimits', maplim);
    title("tRNS", "FontSize",13)
    nexttile;
    tACS = mean(power_changes.values(tACSselect, :), 1);
    topoplot(tACS, ch32Locations, 'maplimits', maplim);
    title("tACS", "FontSize",13)
    cb = colorbar('Limits', maplim);
    cb.Layout.Tile = 'south';
    cb.FontSize = 11;
    xlabel(t, "Alpha Power Change (%)");
    hold off
end

    

function chan_output = map_chan(chan_input)
% Maps the index or name of a channel to the corresponding name or index of
% the channel, respectively
    global chan_map
    chan_output = zeros(length(chan_input), 1);
    if isa(chan_input, 'string')
        for ch = 1:length(chan_input)
            chan_output(ch) = find(contains(chan_map, chan_input(ch)));
        end
    else
        chan_output = chan_map(chan_input);
    end
end