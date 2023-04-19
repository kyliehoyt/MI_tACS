clear
close all
clc

%% Load Subject Data
Subject_12 = load('Subj012.mat', 'sub').sub;
Subject_16 = load('Subj016.mat', 'sub').sub;
Subject_17 = load('Subj017.mat', 'sub').sub;
fs = Subject_12.Pre.restingState.run.header.SampleRate;
global chan_map
chan_map = string(Subject_12.Pre.restingState.run.header.Label(1:32));

tRNS_Subjects = Subject_12;
tACS_Subjects = [Subject_16, Subject_17];
All_Subjects = [Subject_12, Subject_16, Subject_17];

% Preprocessing
All_Subjects = cleanSubjects(All_Subjects, fs, 0, [1 100]);

%% Compute Performance Metrics
[all_perf,all_to] = compute_performance(All_Subjects, fs);
create_bar_plot(all_perf)

%% Compute Absolute Alpha Power
[ffts, f] = getSubjectFFTs(All_Subjects, fs);

D = fittype('a+b/f.^c','dependent', {'D'}, 'independent', {'f'}, 'coefficients', {'a', 'b', 'c'});
C3 = map_chan("C3");
C4 = map_chan("C4");
[C3_baselines, C3_alpha_power, ~, ~] = fitBaselines(ffts, f, D, C3, [7, 7; 4, 4; 7, 7], [29, 29; 32, 30; 25, 28]);
[C4_baselines, C4_alpha_power, ~, ~] = fitBaselines(ffts, f, D, C4, [8, 8; 3, 5; 7, 7], [28, 30; 32, 32; 25, 25]);

%% Statistical Analysis of BCI Perf Metrics
% Compute statistical significance of CDA and TO pre- and post- stimulation
pre = 1; post = 2;
Pcda = evalBciPerfSignificance(all_perf(:,:,pre), all_perf(:,:,post));
Pto = evalBciPerfSignificance(all_to(:,:,pre), all_to(:,:,post));

% Compute Pearson correlation between CDA and Absolute Alpha Power of C3 and C4 features
[rC3, pC3] = corrAlphaAndCDA(C3_alpha_power, all_perf, "C3");
[rC4, pC4] = corrAlphaAndCDA(C4_alpha_power, all_perf, "C4");

%% Functions for Performance
function [accuracy, timeout_rate] = compute_performance(subjects, fs)
    n_sub = length(subjects);
    sessions = ["Pre", "Post"];
    n_sess = length(sessions);
    n_run = 3;
    [accuracy, timeout_rate] = deal(zeros(n_sub, n_run, n_sess));
    for sub = 1:n_sub
        for sess = 1:n_sess
            session = sessions(sess);
            for run = 1:n_run
                [n_hit, n_miss, n_to] = deal(0);
                trigs = getfield(subjects, {sub}, session, 'Online', 'run', {run}, 'header','EVENT', 'TYP');
                pos = getfield(subjects, {sub}, session, 'Online', 'run', {run}, 'header','EVENT', 'POS');
                task_start = pos(find(trigs==1000)+2);
                trial_outcome = trigs(find(trigs==1000)+3);
                trial_end = pos(find(trigs==1000)+3);
                for t = 1:length(trial_outcome)
                    if trial_outcome(t) == 7693 || trial_outcome(t)==7703
                        n_hit = n_hit+1;
                    elseif trial_end(t)-task_start(t) < 7*fs
                        n_miss = n_miss+1;
                    else
                        n_to = n_to+1;
                    end
                end
                accuracy(sub, run, sess) = n_hit/(n_hit+n_miss);
                timeout_rate(sub, run, sess) = n_to/length(trial_outcome);
            end
        end
    end
end

function create_bar_plot(data)
    avg_data = squeeze(mean(data, 2));
    sem_data = std(data, 0, 2)/sqrt(size(data, 2));
    xlab = categorical({'Pre', 'Post'});
    xlab = reordercats(xlab,{'Pre', 'Post'});
    h=bar(xlab, avg_data');
    hold on
    errorbar(h(1).XEndPoints,avg_data(1, :),sem_data(1, :),'LineStyle','none','Color','k','LineWidth',2)
    errorbar(h(2).XEndPoints,avg_data(2, :),sem_data(2, :),'LineStyle','none','Color','k','LineWidth',2)
    errorbar(h(3).XEndPoints,avg_data(3, :),sem_data(3, :),'LineStyle','none','Color','k','LineWidth',2)
    for sub = 1:3
        for s = 1:2
            scatter(repmat(h(sub).XEndPoints(s), 3, 1), squeeze(data(sub,:, s)) ...
            ,10,'MarkerFaceColor', h(sub).FaceColor, 'MarkerEdgeColor','k','LineWidth' ...
            ,1,'XJitter','randn','XJitterWidth',.05)
        end
    end
    leg = legend(h, ["Sub 12 - tRNS" "Sub 16 - tACS" "Sub 17 - tACS"]);
    leg.FontSize = 11;
    ylabel("Command Delivery Accuracy", 'FontSize',13)
    xlabel("Session", "FontSize",13)
    ax = gca;
    ax.FontSize = 11;
end

function plotFFT_dB(f, fft_dB)
    plot(f,fft_dB) 
    title("Single-Sided Power Spectrum")
    xlabel("f (Hz)")
    ylabel("Power (dB)")
    xlim([0 50])
end


function [subject_FFT, f] = getSubjectFFTs(subjects, fs)
    n_sub = length(subjects);
    sessions = {'Pre', 'Post'};
    n_sess = length(sessions);
    c = cell(n_sess, n_sub);
    subject_FFT = cell2struct(c, sessions);
    for sub = 1:n_sub
        for sess = 1:n_sess
            session = sessions{sess};
            period = getfield(subjects, {sub}, session, 'restingState', 'run', 'eeg');
            [pxx, f] = pwelch(period(10*fs:end-10*fs,:), fs, fs/2, 1:50, fs); % middle 40 seconds of resting period
            subject_FFT = setfield(subject_FFT, {sub}, session, 10*log(pxx)); % 1 second segments
        end
    end
end


% function to compute the absolute alpha power: look at paper, see how to
% take out 1/f noise to subtract the unknown baseline
function [baselines, absolute_alpha_power, lower_a, upper_b] = fitBaselines(ffts, f, model, ch_num, lower_a, upper_b)
    n_sub = length(ffts);
    sessions = {'Pre', 'Post'};
    n_sess = length(sessions);
    c = cell(n_sess, n_sub);
    baselines = cell2struct(c, sessions);
    absolute_alpha_power = zeros(n_sub, n_sess);
    fo = fitoptions('StartPoint', [-30, 150, 1], 'Method', 'NonlinearLeastSquares');
    subject_names = ["12", "16", "17"];
    for sub = 1:n_sub
        figure();
        h = zeros(2,1);
        for sess = 1:n_sess
            hold off
            session = sessions{sess};
            full_fft = getfield(ffts, {sub}, session);
            ch_fft = full_fft(:, ch_num); % only 1 channel
            % get bounds
            if nargin<=4
                psd_fig = figure();
                plotFFT_dB(f, ch_fft)
                lower_a(sub, sess) = input("Enter the lower frequency of the alpha peak: ");
                upper_b(sub, sess) = input("Enter the upper frequency of the beta peak: ");
                close(psd_fig);
            end
            % apply bounds
            fitting_fft = vertcat(ch_fft(1:lower_a(sub, sess)), ch_fft(upper_b(sub, sess):end));
            fitting_f = vertcat(f(1:lower_a(sub, sess)), f(upper_b(sub, sess):end));
            % fit curve to get baseline
            fitobject = fit(fitting_f, fitting_fft, model, fo);
            baselines = setfield(baselines, {sub}, session, fitobject);
            % absolute alpha power
            alpha_fft = ch_fft(lower_a(sub, sess):13);
            [alpha_peak, idx] = max(alpha_fft);
            baseline_at_peak = fitobject(idx+lower_a(sub, sess)-1);
            absolute_alpha_power(sub, sess) = alpha_peak - baseline_at_peak;
            % plot
            h(sess) = subplot(1, 2, sess);
            ax = gca;
            lines = plot(fitobject,f, ch_fft);
            set(lines(2), 'LineWidth', 1.5);
            hold on
            xlin(idx+lower_a(sub, sess)-1, string(absolute_alpha_power(sub, sess)), baseline_at_peak, alpha_peak, alpha_peak);
            title(ax, session + " Stimulation", 'FontSize', 13)
            xlabel(ax, "Frequency (Hz)", "FontSize", 13)
            ylabel(ax, "Power (dB)", "FontSize", 13)
            xlim([0 50])
            ax.FontSize = 11;
            leg = legend({'PSD', 'Baseline'});
            leg.FontSize = 11;
        end
        linkaxes(h, 'xy');
        ylim('auto')
        sgtitle(map_chan(ch_num) + " PSD and Baseline Noise of Subject " + subject_names(sub), 'FontSize', 14)
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

function [hl,ht] = xlin(x,txt,ylo,yhi,ytxt) 
% Documentation: 
%   x       = x-Position
%   txt     = Text String
%   ylo     = Low y-Value (Start)
%   yhi     = High y-Value (End)
%   ytxt    = Text Starting Position
hl = plot([x x],[ylo yhi],'DisplayName',txt, 'LineWidth',1);
ht = text(x,ytxt, txt, 'Horiz','left', 'Vert','bottom');
end


% function to see if there is a statistically significant difference in CDA and TO before and after
% tACS
function P = evalBciPerfSignificance(pre, post)
    P = zeros(2,1);
    tRNSPre = pre(1,:);
    tACSPre = pre(2:3,:);
    tRNSPost = post(1,:);
    tACSPost = post(2:3,:);

    [~, P(1)] = ttest2(tRNSPre, tRNSPost); 
    [~, P(2)] = ttest2(tACSPre(:), tACSPost(:));
end

% function to compute correlation absolute alpha power and CDA
function [r, p] = corrAlphaAndCDA(absolute_alpha_power, accuracy, ch_label)
    n_sub = 3;
    n_sess = 2;
    n_run = 3;
    session_accuracy = zeros(n_sub, n_sess);

    for sub = 1:n_sub
        for sess = 1:n_sess
            session_accuracy(sub, sess) = sum(accuracy(sub,:,sess))/n_run;        
        end
    end

    figure();
    scatter(absolute_alpha_power(:), session_accuracy(:), 40, 'blue', 'filled');
    hold on
    line = plot(fit(absolute_alpha_power(:), session_accuracy(:), 'poly1'));
    line.LineWidth = 2;
    line.LineStyle = '--';
    hold off
    xlabel("Absolute Alpha Power (dB)", 'Fontsize', 13);
    ylabel("Command Delivery Accuracy", 'Fontsize', 13);
    title(sprintf("CDA vs Absolute Alpha Power of %s", ch_label), 'Fontsize', 14);
    ax = gca;
    ax.FontSize = 11;

    [R, P] = corrcoef(absolute_alpha_power, session_accuracy);
    r = R(1, 2);
    p = P(1, 2);
end

