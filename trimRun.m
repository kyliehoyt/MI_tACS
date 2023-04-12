function trimmed_run = trimRun(run, n_chan)
    last_row = find(run(:, 1), 1, "last");
    trimmed_run = run(1:last_row, 1:n_chan);
end