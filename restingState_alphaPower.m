clear
close all
clc
%% Load subject data
Subject_12 = load('Subj012.mat', 'sub').sub;
Subject_16 = load('Subj016.mat', 'sub').sub;
Subject_17 = load('Subj017.mat', 'sub').sub;
fs = Subject_12.Pre.restingState.run.header.SampleRate;


% Subject groups
% tRNS_Subjects = Subject_12;
% tACS_Subjects = [Subject_16, Subject_17];
All_Subjects = [Subject_12, Subject_16, Subject_17];

% Global variables
global chan_map
chan_map = string(Subject_12.Pre.restingState.run.header.Label(1:32));
win = 1; % seconds
winfreq = 16; % Hz
alpha_beta = [4 30]; % Hz
alpha = [8 13]; % Hz

%% fft and anova method
rs_alpha = buildPowerTable(All_Subjects, fs, ["Sham", "tACS", "tACS"]);
norm_rs_alpha = normalizeAlphaPower(rs_alpha);
[aov_results, signif_chans] = rsPowerAOV(norm_rs_alpha);
plot_chan = map_chan("C4");
createPowerBarPlot(norm_rs_alpha, plot_chan);
% Some_Subjects aov: 1, 2, 12, 15, 21, 27
% All_Subjects aov: 15, 21
%% PSD and t-test method (no longer using)
% tRNS_rest_psd = restruns2psd(tRNS_Subjects, fs, alpha_beta, 2, win, winfreq);
% tACS_rest_psd = restruns2psd(tACS_Subjects, fs, alpha_beta, 2, win, winfreq);
% All_rest_psd = restruns2psd(All_Subjects, fs, alpha_beta, 2, win, winfreq);
% tRNS_rest_psd_pvals = psd_ttest(tRNS_rest_psd, alpha, true);
% tACS_rest_psd_pvals = psd_ttest(tACS_rest_psd, alpha, false);
% All_rest_psd_pvals = psd_ttest(All_rest_psd, alpha, true);

%% Method using PSD and t-test (No longer using)
function rest_psd = restruns2psd(subjects, fs, psd_band, resolution, win, winfreq)
    n_sub = length(subjects);
    sessions = ["Pre", "Post"];
    n_sess = length(sessions);
    n_chan = 32;
    psd_min = psd_band(1)/resolution + 1; % 4 Hz to bin 3
    psd_max = psd_band(2)/resolution + 1; % 30 Hz to bin 16
    rest_psd = zeros(n_sub, n_sess, psd_max-psd_min+1, n_chan);
    WSize = floor(win*fs);	    % length of each data frame
    nOlap = floor(WSize/winfreq);  % overlap of successive frames
    for sub = 1:n_sub
        for sess = 1:n_sess
            session = sessions(sess);
            run = getfield(subjects, {sub}, session, 'restingState', 'run', 'eeg');
            run = run(:, 1:n_chan);
            [p,~] = pwelch(run, WSize, nOlap, fs/resolution, fs);
            rest_psd(sub, sess, :, :) = p(psd_min:psd_max, :);
        end
    end
end

function p_val = psd_ttest(psds, test_band, grouped)
    %Paired t-test between pre and post resting state power in the test
    %band. Each feature from pre and post PSD (4-30 Hz x 32 channels) form a pair.
    n_sub = size(psds, 1);
    band_min = floor(test_band(1)/2) + 1; % 8 Hz to bin 5
    band_max = ceil(test_band(2)/2) + 1; % 12 Hz to bin 7
    n_feat = (band_max-band_min+1)*size(psds, 4);
    if grouped
        pre_psd = reshape(squeeze(psds(:, 1, band_min:band_max, :)), 1, n_sub*n_feat);
        post_psd = reshape(squeeze(psds(:, 2, band_min:band_max, :)), 1, n_sub*n_feat);
        [~, p_val] = ttest(pre_psd, post_psd);
    else
        p_val = zeros(n_sub, 1); 
        for sub = 1:n_sub
            pre_psd = reshape(squeeze(psds(sub, 1, band_min:band_max, :)), 1, n_feat);
            post_psd = reshape(squeeze(psds(sub, 2, band_min:band_max, :)), 1, n_feat);
            [~, p_val(sub)] = ttest(pre_psd, post_psd);
        end
    end
end

%% Method using fft and two-way anova (from paper)
function mean_alpha_power = meanAlphaPower(segment, fs, f_band)
    L = size(segment, 1);
%     freq_dom = fft(segment, L, 1);
%     two_sided = abs(freq_dom/L).^2;
%     single_sided = two_sided(1:L/2+1, :);
%     f = fs*(0:(L/2))/L;
    %[single_sided, f] = pwelch(period(10*fs:end-10*fs,:), fs, fs/2, 1:50, fs);
    single_sided_alpha = single_sided(find(f==f_band(1),1):find(f==f_band(2),1), :);
    mean_alpha_power = mean(single_sided_alpha, 1);
end

function segment_means = getSegmentPowers(period, fs, seglen)
    L = size(period, 1);
    wSize = floor(seglen*fs); % seconds to samples
    n_seg = floor(L/wSize);
    segment_means = zeros(n_seg, size(period, 2));
    for seg = 1:n_seg
        segment = period((seg-1)*wSize+1:(seg)*wSize, :);
        segment_means(seg, :) = meanAlphaPower(segment, fs, [8 13]);
    end
end

function rs_alpha_power = buildPowerTable(subjects, fs, treatment_labels)
    n_sub = length(subjects);
    sessions = ["Pre", "Post"];
    n_sess = length(sessions);
    n_chan = 32;
    [time, treatment, segment_means] = deal([]);
    for sub = 1:n_sub
        for sess = 1:n_sess
            session = sessions(sess);
            period = getfield(subjects, {sub}, session, 'restingState', 'run', 'eeg');
            period = cleanPeriod(period, n_chan);
            new_means = getSegmentPowers(period, fs, 1); % 1 second segments
            segment_means = vertcat(segment_means, new_means);
            treatment = vertcat(treatment, repmat([treatment_labels(sub)], size(new_means, 1), 1));
            time = vertcat(time, repmat(session, size(new_means, 1), 1));
        end
    end
    rs_alpha_power = struct('values', segment_means, 'treatment', treatment, 'time', time);
end

function [aov_results, signif_chans] = rsPowerAOV(power_struct)
    n_chan = 32;
    fields = {'p_treatment', 'p_time', 'p_interaction', 'tbl', 'stats'};
    c = cell(length(fields), n_chan);
    aov_results = cell2struct(c, fields);
    for chan = 1:n_chan
        [p, tbl, stats] = anovan(power_struct.values(:, chan), {power_struct.treatment, power_struct.time}, 'model', 2, 'varnames', {'treatment', 'time'});
        aov_results(chan).chan_name = map_chan(chan);
        aov_results(chan).p_treatment = p(1);
        aov_results(chan).p_time = p(2);
        aov_results(chan).p_interaction = p(3);
        aov_results(chan).tbl = tbl;
        aov_results(chan).stats = stats;
    end
    close all hidden
    signif_chans.chan = map_chan([aov_results.p_interaction] < 0.05);
    signif_chans.chan_num = find([aov_results.p_interaction] < 0.05);
end

function norm_power_struct = normalizeAlphaPower(power_struct)
    pre_sham_powers = power_struct.values(find(power_struct.treatment == "Sham" & power_struct.time == "Pre"), :);
    pre_sham_mean = mean(pre_sham_powers, 1);
    power_struct.values = power_struct.values./pre_sham_mean.*100;
    norm_power_struct = power_struct;
end

function clean_period = cleanPeriod(period, n_chan)
    last_row = find(period(:, 1), 1, "last");
    clean_period = period(1:last_row, 1:n_chan);
end

function createPowerBarPlot(power_struct, chan_num)
    treatments = ["Sham", "tACS"];
    n_treat = length(treatments);
    for ch = 1:length(chan_num)
        chan = chan_num(ch);
        figure(ch);
        ax = gobjects(n_treat, 1);
        power_struct.values(find(power_struct.treatment == "Sham" & power_struct.time == "Pre"), :) = deal(100);
        for treat = 1:n_treat
            treatment = treatments(treat);
            ax(treat) = subplot(1, 2, treat);
            pre = squeeze(power_struct.values(find(power_struct.treatment == treatment & power_struct.time == "Pre"), chan));
            post = squeeze(power_struct.values(find(power_struct.treatment == treatment & power_struct.time == "Post"), chan));
            means = [mean(pre), mean(post)];
            sems = [std(pre)/sqrt(length(pre)), std(post)/sqrt(length(post))];
            xlab = categorical({'Pre', 'Post'});
            xlab = reordercats(xlab,{'Pre', 'Post'});
            r = bar(xlab(1), means(1), 'FaceColor','r');
            hold on
            o = bar(xlab(2), means(2), 'FaceColor','g');
            errorbar(r.XEndPoints,means(1), sems(1),'LineStyle','none','Color','k','LineWidth',2)
            errorbar(o.XEndPoints,means(2), sems(2),'LineStyle','none','Color','k','LineWidth',2)
            xlabel(treatment);
            ylabel("Alpha Band Power (%)")
            axis(ax(treat), 'tight')
        end
        sgtitle("Relative Alpha Power at " + map_chan(chan))
        linkaxes(ax, 'y')
        hold off
    end
end
       
function chan_output = map_chan(chan_input)
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