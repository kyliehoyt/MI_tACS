clear
close all
clc

Subject_12 = load('Subj012.mat', 'sub').sub;
Subject_16 = load('Subj016.mat', 'sub').sub;
Subject_17 = load('Subj017.mat', 'sub').sub;
fs = Subject_12.Pre.restingState.run.header.SampleRate;

tRNS_Subjects = Subject_12;
tACS_Subjects = [Subject_16, Subject_17];
All_Subjects = [Subject_12, Subject_16, Subject_17];

%strrep(inputname(1), '_', ' ')
win = 1; % seconds
winfreq = 16; % Hz
alpha_beta = [4 30];
alpha = [8 13];

tRNS_rest_psd = restruns2psd(tRNS_Subjects, fs, alpha_beta, 2, win, winfreq);
tACS_rest_psd = restruns2psd(tACS_Subjects, fs, alpha_beta, 2, win, winfreq);
All_rest_psd = restruns2psd(All_Subjects, fs, alpha_beta, 2, win, winfreq);
tRNS_rest_psd_pvals = psd_ttest(tRNS_rest_psd, alpha, true);
tACS_rest_psd_pvals = psd_ttest(tACS_rest_psd, alpha, false);
All_rest_psd_pvals = psd_ttest(All_rest_psd, alpha, true);

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
    % Paired t-test between pre and post resting state power in the test
    % band. Each feature from pre and post PSD (4-30 Hz x 32 channels) form a pair.
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
