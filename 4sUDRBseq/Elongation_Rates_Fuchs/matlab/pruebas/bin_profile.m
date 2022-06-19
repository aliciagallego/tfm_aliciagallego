function res = bin_profile (profile, bin_size)
    number_of_bins = ceil(length(profile)/bin_size);
    length_padding = mod(bin_size-mod(length(profile),bin_size),bin_size);
    res = nanmean(reshape([profile' nan(1,length_padding)],bin_size,number_of_bins));
end