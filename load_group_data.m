clc
clear
addpath('C:\Program Files\MATLAB\R2022b\toolbox\')
folderpath = pwd + "\Group_data_2";
addpath(genpath(folderpath));

treatmentgroups = ["tACS", "tRNS"];
subjects = {["016", "017"], "012"};
session_times = ["Pre", "Post"];
recording_types = ["Online", "restingState"];

for group_id = 1:length(treatmentgroups) % ["tACS", "tRNS"]
    for subject = subjects{1, group_id} % {["016", "017"], "012"}
        treatment = treatmentgroups(group_id);
        sub = struct();
        for session_time = session_times % ["Pre", "Post"]
            for session_type = recording_types % ["Online", "restingState"]

                folder = "Group_data_2/" + treatment + "/" + session_time + "/Subject_" + subject + "_tACSMI/Subject_" + subject + "_tACSMI_Session_001/tACSMI_" + session_type + "_Visual";
                folderpath = pwd + "/" + folder;
                addpath(folderpath);

                fileList = dir(fullfile(folder, '*.gdf'));
                for r = 1:length(fileList)
                    [s, h] = sload(fileList(r).name);
                    sub = setfield(sub, session_time, session_type, 'run', {r}, 'eeg', s);
                    sub = setfield(sub, session_time, session_type, 'run', {r}, 'header', h);
                end
            end
        end
        save(strcat("Subj", subject, ".mat"), 'sub');
    end
end

subj12 = load('Subj012.mat', 'sub').sub;
subj16 = load('Subj016.mat', 'sub').sub;
subj17 = load('Subj017.mat', 'sub').sub;

% # Deland's File Structure
% # > <treatmentgroup>
% #     > Post
% #        > Subject_<id>_tACSMI
% #            > Subject_<id>_tACSMI_Session_001
% #                 > tACS_MI_Offline_Visual
% #                 > tACS_MI_Online_Visual
% #                     > 3 .gdf files
% #                 > tACS_MI_restingState_Visual
% #                     > 1 .gdf file
% #     > Pre
% #        > Subject_<id>_tACSMI
% #            > Subject_<id>_tACSMI_Session_001
% #                 > tACS_MI_Offline_Visual
% #                 > tACS_MI_Online_Visual
% #                     > 3 .gdf files
% #                 > tACS_MI_restingState_Visual
% #                     > 1 .gdf file

% # Our structs
% # > subj<id>
% #     > Post
% #        > Online
% #             > run(x)
% #                 > eeg
% #                 > header
% #        > restingState
% #             > run
% #                 > eeg
% #                 > header
% #     > Pre
% #        > Online
% #             > run(x)
% #                 > eeg
% #                 > header
% #        > restingState
% #             > run
% #                 > eeg
% #                 > header