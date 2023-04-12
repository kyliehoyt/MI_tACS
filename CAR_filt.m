function CAR_filtered = CAR_filt(sig)
    dim = min(size(sig));
    Filter = -1/(dim)*ones(dim) + eye(dim);
    CAR_filtered = sig*Filter;
end