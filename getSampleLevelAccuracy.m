function sample_level_accs = getSampleLevelAccuracy()
    addpath('C:\Program Files\MATLAB\R2022b\toolbox\')
    folderpath = pwd + "\Group_data_2";
    addpath(genpath(folderpath));
    
    treatmentgroups = ["tRNS", "tACS"];
    subjects = {"012", ["016", "017"]};
    session_times = ["Pre", "Post"];
    sample_level_accs = zeros(3, 3, 2); % 3 subject x 3 runs x 2 sessions
    alpha = 0.95;
    n_sub = 0;
    for group_id = 1:length(treatmentgroups) % ["tACS", "tRNS"]
        for subject = subjects{1, group_id} % {["016", "017"], "012"}
            treatment = treatmentgroups(group_id);
            n_sub = n_sub + 1;
            for n_sess = 1:length(session_times) % ["Pre", "Post"]
                session_time = session_times(n_sess);
        
                folder = "Group_data_2/" + treatment + "/" + session_time + "/Subject_" + subject + "_tACSMI/Subject_" + subject + "_tACSMI_Session_001/tACSMI_Online_Visual";
                folderpath = pwd + "/" + folder;
                addpath(folderpath);
        
                fileList = dir(fullfile(folder, '*.txt'));
                for r = 1:length(fileList)
                    fid = fileList(r).name;
                    data = importdata(fid);
                    starts = find(data.data(:, 1) == 7701 | data.data(:, 1) == 7691); % left or right
                    starts(end+1) = size(data.data, 1)-1;
                    sample_level_acc_sum = 0;
                    for s = 1:length(starts)-1
                        if data.data(starts(s), 1) == 7691
                            class = 2; % right
                        else
                            class = 3; % left
                        end
                        Postprob = data.data(starts(s)+1:starts(s+1)-2, class);
                        Prior = [50; Postprob(1:end-1)];
                        Prob = (Postprob - alpha.*Prior)./(1-alpha);
                        sample_level_acc_sum = sample_level_acc_sum + sum(Prob>50)/length(Prob);
                    end
                    sample_level_accs(n_sub, r, n_sess) = sample_level_acc_sum/20;
                end
            end
        end
    end
end