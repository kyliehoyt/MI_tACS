function subjects = cleanSubjects(subjects, fs, filter_order, filter_band)
    n_sub = length(subjects);
    session_times = ["Pre", "Post"];
    session_types = ["restingState", "Online"];
    n_time = length(session_times);
    n_type = length(session_types);
    n_chan = 32;
    for i_sub = 1:n_sub
        for i_time = 1:n_time
            session_time = session_times(i_time);
            for i_type = 1:n_type
                session_type = session_types(i_type);
                n_run = length(getfield(subjects, {i_sub}, session_time, session_type, 'run'));
                for i_run = 1:n_run
                    run = getfield(subjects, {i_sub}, session_time, session_type, 'run', {i_run}, 'eeg');
                    trimmed_run = trimRun(run, n_chan);
                    time_trimmed_run = butter_filt(trimmed_run, fs, filter_order, [filter_band(1) filter_band(2)]);
                    space_time_trimmed_run = CAR_filt(time_trimmed_run);
                    subjects = setfield(subjects, {i_sub}, session_time, session_type, 'run', {i_run}, 'eeg', space_time_trimmed_run);
                end
            end
        end
    end
end
