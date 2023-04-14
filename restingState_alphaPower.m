clear
close all
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

%% fft and anova method
norm_rs_alpha = buildPowerTable(All_Subjects, fs, {'tRNS', 'tACS', 'tACS'});
[aov_results, signif_chans] = rsPowerAOV(norm_rs_alpha);
plot_chan = map_chan(["C3", "C4"]);
createPowerBarPlot(norm_rs_alpha, plot_chan);
createTopoplots(norm_rs_alpha, ch32Locations);


%% Method using fft and two-way anova (from paper)

function rs_alpha_power = buildPowerTable(subjects, fs, treatment_labels)
    n_sub = length(subjects);
    sessions = {'Pre', 'Post'};
    n_sess = length(sessions);
    n_chan = 32;
    [time, treatment] = deal(cell(n_sub*n_sess, 1));
    segment_means = zeros(n_sub*n_sess, n_chan);
    for sub = 1:n_sub
        for sess = 1:n_sess
            session = sessions{sess};
            period = getfield(subjects, {sub}, session, 'restingState', 'run', 'eeg');
            period = trimRun(period, n_chan);
            [pxx, ~] = pwelch(period(10*fs:end-10*fs,:), fs, fs/2, 1:50, fs);
            mean_alpha_powers = sum(pxx(8:12, :), 1);
            element = (sub-1)*n_sess + sess;
            segment_means(element, :) = mean_alpha_powers;
            treatment{element, 1} = treatment_labels{sub};
            time{element, 1} = session;
        end
    end
    rs_alpha_power = struct('values', segment_means, 'treatment', {treatment}, 'time', {time});
    rs_alpha_power = normalizeAlphaPower(rs_alpha_power);
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

function power_struct = normalizeAlphaPower(power_struct)
    Pre = cellfun(@(m)isequal(m,"Pre"),power_struct.time);
    Post = cellfun(@(m)isequal(m,"Post"),power_struct.time);
    pre_powers = power_struct.values(Pre, :);
    power_struct.values(Post,:) = power_struct.values(Post,:)./pre_powers.*100;
    power_struct.values(Pre,:) = power_struct.values(Pre,:)./pre_powers.*100;
end


function createPowerBarPlot(power_struct, chan_num)
    treatments = ["tRNS", "tACS"];
    n_treat = length(treatments);
    for ch = 1:length(chan_num)
        chan = chan_num(ch);
        figure(ch);
        ax = gobjects(n_treat, 1);
        for treat = 1:n_treat
            treatment = treatments(treat);
            ax(treat) = subplot(1, 2, treat);
            treatmentselect = cellfun(@(m)isequal(m,treatment),power_struct.treatment);
            preselect = cellfun(@(m)isequal(m,"Pre"),power_struct.time);
            postselect = cellfun(@(m)isequal(m,"Post"),power_struct.time);
            pre = squeeze(power_struct.values(treatmentselect & preselect, chan));
            post = squeeze(power_struct.values(treatmentselect & postselect, chan));
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

function createTopoplots(power_struct, ch32Locations)
    tRNSselect = cellfun(@(m)isequal(m,"tRNS"),power_struct.treatment);
    tACSselect = cellfun(@(m)isequal(m,"tACS"),power_struct.treatment);
    preselect = cellfun(@(m)isequal(m,"Pre"),power_struct.time);
    postselect = cellfun(@(m)isequal(m,"Post"),power_struct.time);
    maplim = [50 200];
    figure();
    t = tiledlayout(1,4, 'TileSpacing','compact');
    nexttile;
    pre_tRNS = mean(power_struct.values(tRNSselect & preselect, :), 1);
    topoplot(pre_tRNS, ch32Locations, 'maplimits', maplim);
    title("Pre tRNS", "FontSize",13)
    nexttile;
    post_tRNS = mean(power_struct.values(tRNSselect & postselect, :), 1);
    topoplot(post_tRNS, ch32Locations, 'maplimits', maplim);
    title("Post tRNS", "FontSize",13)
    nexttile;
    pre_tACS = mean(power_struct.values(tACSselect & preselect, :), 1);
    topoplot(pre_tACS, ch32Locations, 'maplimits', maplim);
    title("Pre tACS", "FontSize",13)
    nexttile;
    post_tACS = mean(power_struct.values(tACSselect & postselect, :), 1);
    topoplot(post_tACS, ch32Locations, 'maplimits', maplim);
    title("Post tACS", "FontSize",13)
    cb = colorbar('Limits', maplim);
    cb.Layout.Tile = 'south';
    xlabel(t, "Power Change (%)");
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