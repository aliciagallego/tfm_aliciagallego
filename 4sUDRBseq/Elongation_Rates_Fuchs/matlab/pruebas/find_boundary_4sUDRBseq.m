function ind = find_boundary_4suDRBseq(tss_pos, exons_pos, chr_data, norm_param, bin_size)
% This function detect the boundary of a 4sUDRB-seq experiment, it gets as
% input the signal of the experiments for a DRB treated cells and of cells
% after a few minutes of DRB removal and return the predicted boundary
% Input variables:
% tss_pos - array of 4 integers
%   1.chromosome (index to the chr_data)
%   2.strand: 1 - positive 2 - negative
%   3.smaller edge of the transcript (index-wise not the 5' end)
%   4.larger edge of the transcript (index-wise not the 3' end)
% For example: [ 18           1     9102627     9134343]
% exons_pos - a cell array of arrays of exons positions
% For example: {[9102627 9102795],	[9117835 9117901],	[9119323 9119386]};
% chr_data - cell array with size (number of chromosomes) x 2
% in chr_data{c,i} we will have the binned data for chromosme index c for
% experiment i, where i = 1 is the zero (DRB treated experiment) and i = 2
% is data after DRB removal. The data for each chromosome is an array of
% the chromosome size divided by the bin size.
% norm_param - normalization parameter, this is due to having different
% amount of reads for the two experiments, so norm_param will be a two
% element array, in place i (1,2) it will have the sum of all signal in the
% chromosomes of chr_data for experiment i.
% bin_size - by which the chr_data was built

% Extracting the intron data of this gene
data_gene_introns = ...
    extract_intron_data(chr_data, tss_pos, exons_pos, bin_size);

% Checking that we have enough data points
if min(sum(~isnan(data_gene_introns) & (data_gene_introns > 0),2)) > 20
    
    % Normalizing by the experiment total signal
    data = data_gene_introns ./ repmat(norm_param', 1, size(data_gene_introns,2));
    raw_data = data;
    
    
    % As we use the 0-time point for normalization we will add
    % pseudocounts to it (so we won't divide by zero)
    data(1,:) = data(1,:) + 1/(bin_size*norm_param(1));

    % Smoothing the data with cubic spline
    data = smooth_data(data);
    
    % If a gene is longer than 60K we take only the begining
    % Notice that this is specificly tuned to 4/8 minutes. Longer times
    % might need more than 60Kb.
    data = data(:,1:min(size(data,2), ceil(60000 / bin_size)));
    disp(size(data));
  

   % plot(data);
    % For each point we will sum 5Kb (with weights) upstream it. As the
    % phenomenoon we are looking for is an "area" one and not a specific
    % point one
    % ALICIA_PLAY WITH DIFFERENT VALUES THAN 5000 
    winsize = ceil(5000 / bin_size);
    data_summed = cummulative_sum_data(data, winsize);
    
    % As we want to see a relative increase, we would like to divide each
    % point by the "abundent" value from this point on. This is a "trick"
    % that will give us an estimate on the "background" level of signal in
    % this specific gene.
    data_div = data_divide_by_histpeak(data, data_summed, winsize);
    
    % We will now find the first aproximation: where the data_div of the
    % wave data is lower than the '0' profile by some threshold (the most
    % upstream location the fulfills this crateria)
    threshold = 1.1; % 10% close to the 0 time-point
    ind = nan;
    
    x = (data_div(2,:)+0.01)./(data_div(1,:)+0.01);
    x(isinf(x)) = nan;

    
    if i ~= 1
        % The "minimum" upstream location for peak to be found is "winsize"
        % as we include data for winsize points upstream in each point,
        % we are not senstitive enough to locations upstream winsize.
        cur_ind = find(x(winsize:end)<threshold,1)+winsize-1;
        if ~isempty(cur_ind)
            % as we sum winsize data upstream the value, we take half
            % win_size in our estimate
            ind = cur_ind - winsize/2;
        end
    end
    
    % optional filtering
    if ind <= winsize/2
        ind = nan;
    end
    
    % As the previous method finds a rough estimate of the boundary point
    % we will try to find a more subtle breaking point in a different
    % more delicate procedure
    ind = subtle_estimate(ind, raw_data, bin_size);
else
    ind = nan;
end
ind = ind * bin_size;
end


function val = histpeak(data, hist_bins)
% Find the peak of the histogram (the most common value)
% hist_bins - number of bins for the histogram
[dist bins_values] = hist(data, hist_bins);
val = nanmean(bins_values(dist == max(dist)));
end

function data_summed = cummulative_sum_data(data, winsize)
filter2use = fspecial('gaussian',winsize*2,winsize/2);
filter2use = filter2use(winsize,:);
filter2use((winsize+1):end) = 0;
filter2use = filter2use / sum(filter2use);

data_no_nan = data;
data_no_nan(isnan(data_no_nan)) = 0;
data_summed = filter2(filter2use,data_no_nan) ./ filter2(filter2use, ~isnan(data));
end

function data_div = data_divide_by_histpeak(data, data_summed, winsize)
data_div = nan(2,size(data,2));
peaks = nan(2,length(ceil(winsize/2):size(data,2)));

% % % Running on all the location from a certain point on
index_in_peaks = 0;
for i=ceil(winsize/2):size(data,2)
    index_in_peaks = index_in_peaks + 1;
    for j=1:2
        cur_data = data(j,(i+1):end);
        cur_data = cur_data(cur_data > 0);
        % We caluculate the hist peak only if we have at least 10 data
        % points in the relevant area
        if length(cur_data) < 10
            % If we don't have enough data we take the last peak or the
            % average of signal
            if index_in_peaks > 1
                peaks(j,index_in_peaks) = peaks(j,index_in_peaks-1);
            else
                peaks(j,index_in_peaks) = nanmean(data(j,:))+0.01;
            end
        else
            peaks(j,index_in_peaks) = histpeak(cur_data,15);
        end
        data_div(j,i) = data_summed(j,i) ./ peaks(j,index_in_peaks);
    end
end
end

function data = smooth_data(data)
pc = 0.00001;
for i=1:2
    cur_indices = find(~isnan(data(i,:)));
    data(i,:) = csaps(cur_indices,data(i,cur_indices),pc,1:size(data,2));
end
end


function new_ind  = subtle_estimate(ind, data, bin_size)
% After having an initial estimation of the bounderies we will try to
% improve the estimation using information from the pattern of the wave.
% More specifically we will try to see when the derivative will be most
% extreme (so we will know it is decreasing) and then find the place where
% it gets to a plateau and this will be the boundary.

% Here we are going to use the peakfinder built by N. Yoder (mathworks)
% http://www.mathworks.com/matlabcentral/fileexchange/25500-peakfinder

if ~isnan(ind)
    pc = 0.00001; % Parameter for the cubic spline
    for i=1:2
        cur_indices = find(~isnan(data(i,:)));
        % smoothing the data, we will add the '0' profile pseudocounts
        if i~=1
            data(i,:) = csaps(cur_indices,data(i,cur_indices),pc,1:size(data,2));
        else
            data(i,:) = csaps(cur_indices,data(i,cur_indices)+0.01,pc,1:size(data,2));
        end
    end
    
    new_ind = nan;
    % Caluclating the derivative
    cd = diff(data(2,:) ./ data(1,:));
    
    % Finding the local minimums in the derivative (upstream the rough boundary)
    [loc loc_str] = peakfinder(-cd(1:ind),nanmax(-cd(1:ind))/6);
    loc_str = loc_str(loc < ind);
    loc = loc(loc < ind);
    loc = loc(loc_str > 0.002*5); % Threshold on the value of the min
    
    if length(loc) >= 1
        loc = loc(end);
        low_derivative = ... Finding the place after the min which is closest to 0
            find(cd(loc:min(ind+(2000/bin_size),length(cd))) >= -0.002,1);
        if ~isempty(low_derivative)
            new_ind = low_derivative + loc-1;
        end
    end
    % If we couldn't find anything better we will take the older rough
    % estimate
    if isnan(new_ind)
        new_ind = ind;
    end
else
    new_ind = nan;
end
end



function data_gene_introns = extract_intron_data(chr_data, tss_pos, exons_pos, bin_size)
strand = tss_pos(2); % +/- strand
c = tss_pos(1); % Chromosome

if strand == 1
    tss_gene = ceil(tss_pos(3)/bin_size);
    tts_gene = ceil(tss_pos(4)/bin_size);
else
    tss_gene = ceil(tss_pos(4)/bin_size);
    tts_gene = ceil(tss_pos(3)/bin_size);
end

start_loc = min(tss_gene, tts_gene);
end_loc = max(tss_gene, tts_gene);

% retrieve the data from the all-genome data
data_gene_introns = nan(2, (end_loc - start_loc)+1);
for i=1:2
    data_gene_introns(i,:) = chr_data{c,i}(start_loc:end_loc);
end

% Put nans in the location of the exons
for ex_index = 1:length(exons_pos)
    ex_data = sort(exons_pos{ex_index});
    ex_data(1) = max(floor(ex_data(1) / bin_size)-start_loc+1,1);
    ex_data(2) = ceil(ex_data(2) / bin_size)-start_loc+1;
    
    data_gene_introns(ex_data(1):ex_data(2)) = nan; 
end

if strand==2 % If the gene is on the reverse strand, reverse the data
    data_gene_introns = data_gene_introns(:,end:-1:1);
end
end