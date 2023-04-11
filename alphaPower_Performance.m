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

[all_perf,all_to] = compute_performance(All_Subjects, fs);

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

% function to see if there is a statistically significant difference in CDA and TO before and after
% tACS

% function to compute the absolute alpha power: look at paper, see how to
% take out 1/f noise to subtract the unknown baseline

% function to compute correlation absolute alpha power and CDA

