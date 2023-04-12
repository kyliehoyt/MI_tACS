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

%% Compute Performance Metrics
[all_perf,all_to] = compute_performance(All_Subjects, fs);

%% Compute Absolute Alpha Power
test_segment = Subject_12.Post.restingState.run.eeg(1:30000, 11);
test_segment = butter_filt(test_segment, fs, 2, [4 30]);
[fft, f] = plotFFT(test_segment, 512);
hold on
D = fittype('A/(1-abs(4*(x-c)/w))*0.5*(1-cos(2*pi()*x/w))+N*x.^-B','dependent', {'D'}, 'independent', {'x'}, 'coefficients', {'N', 'B', 'A', 'c', 'w'});
fo = fitoptions('StartPoint', [0, 0.45, 0.5, 10, 5], 'Method', 'NonlinearLeastSquares', 'Lower', [0, 0.2, 0, 8, 2], 'Upper', [1, 0.8, 1, 12, 6]);
fitobject = fit(f', fft, D, fo);
plot(fitobject, f, fft)



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

function [single_sided, f] = plotFFT(segment, fs)
    L = size(segment, 1);
    n = 2^nextpow2(L);
    %n = fs*2;
    n = L/2;
    freq_dom = fft(segment, n);
    two_sided = abs(freq_dom/n).^2;
    single_sided = two_sided(2:(n/2)+2);
    f = fs*(0:(n/2))/n;
    plot(f,single_sided) 
    title("Single-Sided Power Spectrum of X(t)")
    xlabel("f (Hz)")
    ylabel("|P(f)|^2")
    xlim([0 50])
end

function sig = butter_filt(Raw_Sig, fs, N, band, btype)
    if length(band) > 1 && nargin < 5
        [B, A] = butter(N, [band(1) band(2)].*(2/fs));
        sig = filtfilt(B, A, Raw_Sig);
    else
        [B, A] = butter(N, band*(2/fs), btype);
        sig = filtfilt(B, A, Raw_Sig);
    end
end
% function to see if there is a statistically significant difference in CDA and TO before and after
% tACS

% function to compute the absolute alpha power: look at paper, see how to
% take out 1/f noise to subtract the unknown baseline

% function to compute correlation absolute alpha power and CDA

