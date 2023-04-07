clear
close all
clc

Subject_12 = load('Subj012.mat', 'sub').sub;
Subject_16 = load('Subj016.mat', 'sub').sub;
Subject_17 = load('Subj017.mat', 'sub').sub;

tRNS_subs = Subject_12;
tACS_subs = [Subject_16, Subject_17];
all_subs = [Subject_12, Subject_16, Subject_17];

%strrep(inputname(1), '_', ' ')