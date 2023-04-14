function sig = butter_filt(Raw_Sig, fs, N, band, btype)
    if N == 0
        sig = Raw_Sig;
        return
    end
    if length(band) > 1 && nargin < 5 
        [B, A] = butter(N/2, [band(1) band(2)].*(2/fs));
        sig = filtfilt(B, A, Raw_Sig);
    else
        [B, A] = butter(N, band*(2/fs), btype);
        sig = filtfilt(B, A, Raw_Sig);
    end
end