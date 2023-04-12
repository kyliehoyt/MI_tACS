clear
close all
clc

%% Load Subject Data
Subject_12 = load('Subj012.mat', 'sub').sub;
Subject_16 = load('Subj016.mat', 'sub').sub;
Subject_17 = load('Subj017.mat', 'sub').sub;
fs = Subject_12.Pre.restingState.run.header.SampleRate;

tRNS_Subjects = Subject_12;
tACS_Subjects = [Subject_16, Subject_17];
All_Subjects = [Subject_12, Subject_16, Subject_17];
All_Subjects = cleanSubjects(All_Subjects, fs, 4, [4 30]);

%% Compute Performance Metrics
[all_perf,all_to] = compute_performance(All_Subjects, fs);

%% Compute Absolute Alpha Power
ffts = getSubjectFFTs(All_Subjects, fs);
GA_ffts = calcGAvgFFT(ffts);
test_fft = GA_ffts(1).Pre(:, 15);
n = fs/2;
f = (fs*(0:(n/2))/n)';

D = fittype('A/(1-abs(4*(x-c)/w))*0.5*(1-cos(2*pi()*x/w))+N*x.^-B','dependent', {'D'}, 'independent', {'x'}, 'coefficients', {'N', 'B', 'A', 'c', 'w'});
fo = fitoptions('StartPoint', [0, 0.45, 0.5, 10, 5], 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0.2, 0, 8, 2], 'Upper', [1, 0.8, 1, 12, 6]);
% fitobject = fit(f, test_fft, D, fo);
% plot(fitobject, f, test_fft)
plotFFT(f, test_fft)


%% Functions for Performance
function [accuracy, timeout_rate] = compute_performance(subjects, fs)
    n_sub = length(subjects);
    sessions = ["Pre", "Post"];
    n_sess = length(sessions);
    n_run = 3;
    [accuracy, timeout_rate] = deal(zeros(n_sub, n_sess, n_run));
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
                accuracy(sub, sess, run) = n_hit/(n_hit+n_miss);
                timeout_rate(sub, sess, run) = n_to/length(trial_outcome);
            end
        end
    end
end

function [single_sided, f] = calcSegmentFFT(segment, fs, Display)
    L = size(segment, 1);
    n = 2^nextpow2(L);
    %n = fs*2;
    n = L/2;
    freq_dom = fft(segment, n, 1);
    two_sided = abs(freq_dom/n).^2;
    single_sided = two_sided(2:(n/2)+2, :);
    f = (fs*(0:(n/2))/n)';
    if Display
        plotFFT(f, single_sided);
    end
end

function plotFFT(f, single_sided)
    plot(f,single_sided) 
    title("Single-Sided Power Spectrum")
    xlabel("f (Hz)")
    ylabel("|P(f)|^2")
    xlim([0 50])
end

function segment_FFT = getSegmentFFTs(period, fs, seglen)
    L = size(period, 1);
    wSize = floor(seglen*fs); % seconds to samples
    n_seg = floor(L/wSize);
    segment_FFT = zeros(n_seg, fs/4+1, size(period, 2));
    for seg = 1:n_seg
        segment = period((seg-1)*wSize+1:(seg)*wSize, :);
        [segment_FFT(seg, :, :), ~] = calcSegmentFFT(segment, fs, false);
    end
end

function subject_FFT = getSubjectFFTs(subjects, fs)
    n_sub = length(subjects);
    sessions = {'Pre', 'Post'};
    n_sess = length(sessions);
    n_chan = 32;
    c = cell(n_sess, n_sub);
    subject_FFT = cell2struct(c, sessions);
    for sub = 1:n_sub
        for sess = 1:n_sess
            session = sessions{sess};
            period = getfield(subjects, {sub}, session, 'restingState', 'run', 'eeg');
            period = trimRun(period, n_chan);
            subject_FFT = setfield(subject_FFT, {sub}, session, getSegmentFFTs(period, fs, 1)); % 1 second segments
        end
    end
end

function grand_average_fft = calcGAvgFFT(subject_FFT)
    n_sub = length(subject_FFT);
    sessions = {'Pre', 'Post'};
    n_sess = length(sessions);
    c = cell(n_sess, n_sub);
    grand_average_fft = cell2struct(c, sessions);
    for sub = 1:n_sub
        for sess = 1:n_sess
            session = sessions{sess};
            session_fft = getfield(subject_FFT, {sub}, session);
            grand_average_fft = setfield(grand_average_fft, {sub}, session, squeeze(mean(session_fft, 1)));
        end
    end
end

% function to see if there is a statistically significant difference in CDA and TO before and after
% tACS

% function to compute the absolute alpha power: look at paper, see how to
% take out 1/f noise to subtract the unknown baseline

% function to compute correlation absolute alpha power and CDA

