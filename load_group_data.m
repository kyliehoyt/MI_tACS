addpath('C:\Program Files\MATLAB\R2022b\toolbox\')
folderpath = pwd + "\Group_data_2";
addpath(genpath(folderpath));

treatmentgroups = ["tACS", "tRNS"];
subjects = {["016", "017"], "012"};
session_times = ["Pre", "Post"];
recording_types = ["Online", "restingState"];


for group_id = 1:length(treatmentgroups) % ["tACS", "tRNS"]
    for time_id = 1:length(session_times) % ["Pre", "Post"]
        for subject = subjects{1, group_id} % {["016", "017"], "012"}
            for session_type_id = 1:length(recording_types) % ["Online", "restingState"]
                treatment = treatmentgroups(group_id);
                session_time = session_times(time_id);
                session_type = recording_types(session_type_id);
                
                folder = "Group_data_2/" + treatment + "/" + session_time + "/Subject_" + subject + "_tACSMI/Subject_" + subject + "_tACSMI_Session_001/tACSMI_" + session_type + "_Visual";
                folderpath = pwd + "/" + folder;
                addpath(folderpath);
    
                save_path = "subject_" + subject + "/" + session_time + "/" + session_type;
                if ~exist(pwd + "/" + save_path, 'dir')
                    mkdir(pwd + "/" + save_path);
                end
                fileList = dir(fullfile(folder, '*.gdf'));
                runnum = 1;
                for f = 1:length(fileList)
                    [s, h] = sload(fileList(f).name);
                    sname = save_path + "/s" + string(runnum);
                    hname = save_path + "/h" + string(runnum);
                    save(sname, "s");
                    save(hname, "h");
                    runnum = runnum+1;
                end
            end
        end
    end

end


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

% # Our File Structure
% # > subject_<id>
% #     > Post
% #        > Online
% #             > s1.mat
% #             > s2.mat
% #             > s3.mat
% #             > h1.mat
% #             > h2.mat
% #             > h3.mat
% #        > restingState
% #             > s1.mat
% #             > h1.mat
% #     > Pre
% #        > Online
% #             > s1.mat
% #             > s2.mat
% #             > s3.mat
% #             > h1.mat
% #             > h2.mat
% #             > h3.mat
% #        > restingState
% #             > s1.mat
% #             > h1.mat