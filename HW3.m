clear all
close all
clc

%% Import Data
subject1 = load('subj1.mat').subj1;
subject2 = load('subj2.mat').subj2;
fs = subject1.offline.run(1).header.fs;
[timeout_rate, accuracies] = deal(zeros(2,6,3));
load("ch32Locations.mat");

%% Compute Performance Metrics
nsub = 1;
for sub = [subject1, subject2]
    for s = 1:6
        for r = 1:3
            n_hit = 0;
            n_miss = 0;
            n_to = 0;
            trigs = sub.online(s).run(r).header.triggers.TYP;
            pos = sub.online(s).run(r).header.triggers.POS;
            task_start = pos(find(trigs==1000)+3);
            trial_out = trigs(find(trigs==1000)+4);
            trial_end = pos(find(trigs==1000)+4);
            for t = 1:20
                if trial_out(t) == 7693 || trial_out(t)==7703
                    n_hit = n_hit+1;
                elseif trial_end(t)-task_start(t) < 7*fs
                    n_miss = n_miss+1;
                else
                    n_to = n_to+1;
                end
            end
            accuracies(nsub, s, r) = n_hit/(n_hit+n_miss);
            timeout_rate(nsub, s, r) = n_to/20;
        end
    end
    nsub = nsub+1;
end

%% Trend in command delivery accuracy over sessions
figure();
[R_acc, p_acc] = create_bar_plot(accuracies);
title("Command Delivery Accuracy")
ylim([0.0, 1.2]);
ylabel("Accuracy")
xticklabels([1, 2, 3, 4, 5, 6])
xlabel("Session")
pairwise_ps_acc = pairwise_ts(accuracies);

%% Trend in timeout rate over sessions
% figure();
% [R_to, p_to] = create_bar_plot(timeout_rate);
% title("Timeout Rate")
% ylabel("Accuracy")
% xticklabels([1, 2, 3, 4, 5, 6])
% xlabel("Session")
% pairwise_ps_to = pairwise_ts(timeout_rate);

%% Feature Discriminability
% Extract features separated by class
[l_feat_1, r_feat_1] = extract_feats(subject1, fs);
[l_feat_2, r_feat_2] = extract_feats(subject2, fs);
% Calculate fisher score for each feature
fishers_s1 = calc_fisher(l_feat_1, r_feat_1);
fishers_s2 = calc_fisher(l_feat_2, r_feat_2);
% Names of features
f_labs = {"4", "6", "8", "10", "12", "14", "16", "18", "20", "22", "24", "26", "28", "30"};
ch_labs = {"Fp1"};
for ch = 2:32
    ch_labs(end+1) = {string(ch32Locations(ch).labels)};
end
% Compute top ten features across sessions
[top_feats, mask] = top_fisher(fishers_s2, f_labs, ch_labs, 10);
top_feats_trend = get_feat_trend(fishers_s2, mask);
% Plot fisher scores over sessions
% plot_fisher(fishers_s2, f_labs, ch_labs, "2")
% fisher_topo(fishers_s2, ch32Locations, "2")
% plot_fisher(fishers_s1, f_labs, ch_labs, "1")
% fisher_topo(fishers_s1, ch32Locations, "1")

%% Trend in top features
% Trends in top feature from session 6 over all sessions
% top feature from session 6
[best_feat_1, best_feat_1_mask] = top_fisher(fishers_s1(6,:,:), f_labs, ch_labs, 1);
[best_feat_2, best_feat_2_mask] = top_fisher(fishers_s2(6,:,:), f_labs, ch_labs, 1);
% trend of feature over sessions
best_feat_trend_1 = get_feat_trend(fishers_s1, best_feat_1_mask);
best_feat_trend_2 = get_feat_trend(fishers_s2, best_feat_2_mask);
% plot trend and compute trendline
feature_trends = vertcat(best_feat_trend_1, best_feat_trend_2);
feature_labels = vertcat(strcat("S1: ", string(best_feat_1)), strcat("S2: ", string(best_feat_2)));
[feat_trend_Rs, feat_trend_ps] = plot_feat_trend(feature_trends, feature_labels);

% Paired t-test between session 1 and session 6 for top feature
pre_post_p = pre_post(feature_trends);

% Correlation of feature discriminability and command delivery accuracy
% over sessions
% top 10 features from session 6
[top10feats_2, top10feats_2_mask] = top_fisher(fishers_s2(6,:,:), f_labs, ch_labs, 10);
% trends of features over sessions
top10feats_trend_2 = get_feat_trend(fishers_s2, top10feats_2_mask);
% session-wise command delivery accuracy
sess_acc_2 = mean(accuracies, 3);
CDA_2 = repmat(squeeze(sess_acc_2(2, :)), 10, 1);
% plot and determine correlation
[corr_rs, corr_ps] = plot_feat_correlation(top10feats_trend_2, CDA_2, top10feats_2);
%% Functions -- Performance Trends
function [Rs, ps] = create_bar_plot(data)
    [Rs, ps] = deal(zeros(3,1)); % (runs, 1)
    xtrend = 1:6;
    avg_data = mean(data, 3);
    std_data = std(data, 0, 3);
    h=bar(avg_data');
    hold on
    h(1).FaceColor='r';
    h(2).FaceColor='b';
    errorbar(h(1).XEndPoints,avg_data(1, :),std_data(1, :),'LineStyle','none','Color','k','LineWidth',2)
    errorbar(h(2).XEndPoints,avg_data(2, :),std_data(2, :),'LineStyle','none','Color','k','LineWidth',2)
    colors = ['r', 'b'];
    for sub = 1:2
        for s = 1:6
            scatter(repmat(h(sub).XEndPoints(s), 3, 1), squeeze(data(sub,s,:)) ...
            ,10,'MarkerFaceColor',colors(sub),'MarkerEdgeColor','k','LineWidth' ...
            ,1,'XJitter','randn','XJitterWidth',.05)
        end
        flatsub = reshape(data(sub, :, :), 1, []);
        lm = polyfitn(repmat(xtrend, 1, 3), flatsub, 1);
        trendline = polyval(lm.Coefficients, xtrend);
        plot(xtrend, trendline, 'Color', colors(sub), 'LineWidth', 2);
        Rs(sub) = lm.R2;
        ps(sub) = lm.p(1);
    end
    allpoints = horzcat(reshape(data(1, :, :), 1, []), reshape(data(2, :, :), 1, []));
    groupedlm = polyfitn(repmat(xtrend, 1, 6), allpoints, 1);
    trendline = polyval(groupedlm.Coefficients, xtrend);
    plot(xtrend, trendline, 'Color', 'k', 'LineWidth',2);
    Rs(3) = groupedlm.R2;
    ps(3) = groupedlm.p(1);
    legend(h, ["Subject 1" "Subject 2"])
end

function [ps] = pairwise_ts(data)
    ps = zeros(2, 5); % (subjects, sessions other than S1)
    for sub = 1:2
        for s = 2:6
            [~, ps(sub,s)] = ttest2(squeeze(data(sub, 1, :)), squeeze(data(sub, s, :)));
        end
    end
end
%% Functions -- Fisher Score 
function [R_feats, L_feats] = extract_feats(subject, fs)
    [R_feats, L_feats] = deal(zeros(6, 30, 14, 32)); % (sessions, trials, freqs, channels)
    WSize = floor(1*fs);	    % length of each data frame
    nOlap = floor(WSize/16);  % overlap of successive frames
    for s = 1:6
        [lt, rt] = deal(1);
        for r = 1:3
            run = subject.online(s).run(r).eeg;
            trigs = subject.online(s).run(r).header.triggers.TYP;
            pos = subject.online(s).run(r).header.triggers.POS;
            class = trigs(find(trigs==1000)+2);
            task_starts = pos(find(trigs==1000)+3);
            task_ends = pos(find(trigs==1000)+4);
            for t = 1:20
                if task_ends(t)-task_starts(t) <= 512
                    [p,~] = pwelch(run(task_starts(t):task_ends(t), :),[], [], [], fs);
                else
                    [p,~] = pwelch(run(task_starts(t):task_ends(t), :),WSize, nOlap, fs/2, fs);
                end
                if class(t) == 769 %Left
                    L_feats(s, lt, :, :) = p(3:16, :);
                    lt = lt+1;
                else
                    R_feats(s, rt, :, :) = p(3:16, :); 
                    rt = rt+1;
                end
            end
        end
    end
end

function fisher_scores = calc_fisher(c1_feats, c2_feats)
    fisher_scores = zeros(6, 14, 32); % (sessions, freqs, channels)
    mean_c1 = squeeze(mean(c1_feats, 2)); % mean across trials in a session
    mean_c2 = squeeze(mean(c2_feats, 2));
    var_c1 = squeeze(var(c1_feats, 0, 2));
    var_c2 = squeeze(var(c2_feats, 0, 2));
    for s = 1:6
        fisher_scores(s, :, :) = abs(mean_c1(s,:,:) - mean_c2(s,:,:))./sqrt(var_c1(s,:,:) + var_c2(s,:,:));
    end
end

function [top_feats, mask] = top_fisher(fishers, f_labs, ch_labs, n_feats)
    rank_weighted_fishers = zeros(14, 32); % (freqs, channels)
    for s = 1:size(fishers,1)
        fish_s = squeeze(fishers(s, :, :));
        [~, idx] = sort(fish_s(:), 'descend');
        rank = zeros(size(fish_s));
        rank(idx) = 1:numel(fish_s);
        rank_weighted_fishers = rank_weighted_fishers + fish_s./rank;
    end
    [~, idx] = sort(rank_weighted_fishers(:), 'descend');
    mask = zeros(size(rank_weighted_fishers));
    mask(idx) = 1:numel(rank_weighted_fishers);
    mask(mask<=n_feats) = 1;
    mask(mask>n_feats) = 0;
    top_feats = {};
    for ch = 1:32
        for f = 1:14
            if mask(f, ch) ~= 0
                top_feats(end+1) = {strcat(string(ch_labs(ch)), ": ", string(f_labs(f)), " Hz")};
            end
        end
    end
end

function plot_fisher(fishers, f_labs, ch_labs, sub_num)
    figure();
    for s = 1:6
        subplot(2, 3, s);
        fisher_sess = squeeze(fishers(s, :, :))';
        fisher_sess = fisher_sess./max(fisher_sess);
        h = heatmap(f_labs, ch_labs, fisher_sess);
        h.XLabel = "Frequency";
        h.YLabel = "Channel";
        h.Title = "Session " + string(s);
        h.ColorLimits = [0 1];
    end
    sgtitle("Fisher Scores -- Subject " + sub_num)
end

function fisher_topo(fishers, ch_locs, sub_num)
    figure();
    for s = 1:6
        subplot(1, 6, s);
        fisher_sess = squeeze(fishers(s, :, :));
        fisher_chan = sum(fisher_sess, 1); % 1x32
        topoplot(fisher_chan, ch_locs,'maplimits', [0.7 3.5]);
        title("Session " + string(s));
        cbar('vert',0, [0.7 3.5]);
    end
    sgtitle("Summed Fisher Score Maps -- Subject " + sub_num);
end
%% Functions - Fisher Score Trends
function feat_trend = get_feat_trend(fishers, mask)
    feat_trend = zeros(sum(mask,"all"), size(fishers,1));
    for s = 1:size(fishers,1)
        masked_feats = squeeze(fishers(s, :,:)).*mask;
        feat_trend(:, s) = masked_feats(masked_feats ~=0)';
    end
end

function [Rs, ps] = plot_feat_trend(trends, feat_labels)
    [Rs, ps] = deal(zeros(size(trends,1), 1)); % (number features, 1)
    xtrend = 1:size(trends, 2); % 1:number of sessions
    figure();
    hold on 
    cmap = hsv(size(trends, 1));
    for f = 1:size(trends, 1)
        scatter(xtrend, squeeze(trends(f, :)),40, cmap(f, :),"filled","o", 'DisplayName', string(feat_labels(f)));
        lm = polyfitn(xtrend, squeeze(trends(f, :)), 1);
        trendline = polyval(lm.Coefficients, xtrend);
        plot(xtrend, trendline, 'LineWidth', 2, 'Color',cmap(f, :), 'DisplayName',strcat(string(feat_labels(f)), " Trend"));
        Rs(f) = lm.R2;
        ps(f) = lm.p(1);
    end
    legend();
    title("Fisher Score Trend over Sessions")
    xlabel("Session")
    xticks(xtrend)
    ylabel("Fisher Score")
end

function p = pre_post(trends)
    [~, p] = ttest(trends(:, 1), trends(:, end));
end

function [rs, ps] = plot_feat_correlation(feat_trends, accs, feat_labels)
    figure();
    hold on
    cmap = hsv(size(feat_trends,1));
    [rs, ps] = deal(zeros(size(feat_trends, 1), 1));
    for f = 1:size(feat_trends, 1)
        scatter(accs(f, :), feat_trends(f, :), 30, cmap(f, :), "filled", "o", 'DisplayName', string(feat_labels(f)))
        [R, P] = corrcoef(accs(f, :), feat_trends(f, :));
        rs(f) = R(1,2);
        ps(f) = P(1,2);
    end
    title("Top 10 Features and Accuracy over Sessions")
    xlabel("Command Delivery Accuracy")
    ylabel("Fisher Score")
    legend('location', 'best');
end
