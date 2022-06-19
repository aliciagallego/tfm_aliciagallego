% cd /media/cc/B/Josemi/TTseq_Feb2022/TTseq_scripts/Elongation_Rates_Fuchs/matlab
% Open writing file and print the header
outfid = fopen('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/Elongation_rate_5min_20220425_20Kb_size_Pull.txt', 'wt');
fprintf(outfid, "Gene_name\tWT_Pull\tTK0_Pull\n");

% WT
all_chr_WT1 = cell(21,2);

% 0 min sample
for i = 1:21
    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/WT0_MG9-11-13/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_WT1{i,1} = bin_profile(raw_data,100);
end

% 5 min sample
for i = 1:21
    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/WT5_MG9-12-14/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_WT1{i,2} = bin_profile(raw_data,100);
end

% H1-TKO
all_chr_TKO1 = cell(21,2);

% 0 min sample
for i = 1:21
    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/TKO0_MG9-15-17/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_TKO1{i,1} = bin_profile(raw_data,100);
end

% 5 min sample
for i = 1:21
    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/TKO5_MG9-16-18/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_TKO1{i,2} = bin_profile(raw_data,100);
end

% For each gene in the Input_genes.txt file
mygenes = tdfread('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_2/4_Input_genes/Input_genes_20Kb.txt', "\t");
tablesize = size(mygenes.Name, 1);

i=1;
while i <= tablesize
    gene_position= [mygenes.Chromosome(i) mygenes.Orientation(i) mygenes.Start(i) mygenes.End(i)];   
    exons=cell(1,mygenes.Exon_number(i));
    exon_starts = strsplit(mygenes.Exon_starts(i,:), ',');
    exon_ends = strsplit(mygenes.Exon_ends(i,:), ',');
    
    j=1;
    while j <= mygenes.Exon_number(i)
        toadd = [str2num(exon_starts{1,j}) str2num(exon_ends{1,j})];
        exons{1,j} = toadd;
        j=j+1;
    end
    
    % The fourth parameter comes from file /media/cc/A/Josemi/PRUEBA_metlab/Elongation_rate/Alignment/RATIOS
    WTbound=find_boundary_4sUDRBseq(gene_position, exons, all_chr_WT1, [1.17 1.35],100);
    TKObound=find_boundary_4sUDRBseq(gene_position, exons, all_chr_TKO1, [1.09 1.00],100);
    
    fprintf(outfid, "%s\t%d\t%d\n", mygenes.Name(i,:), WTbound, TKObound);
    i = i+1;
end

% Close writing file
fclose(outfid);
