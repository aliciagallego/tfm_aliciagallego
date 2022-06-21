% cd /path/Elongation_Rates/matlab
% Open writing file and print the header
outfid = fopen('/path/Elongation_rate_Pull_output.txt', 'wt');
fprintf(outfid, "Gene_name\tWT_Pull\tTK0_Pull\n");

% WT Pull of the replicates
all_chr_WT1 = cell(21,2);

% 0 min
for i = 1:21
    fid = fopen(strcat('/path/WT0_MG9-11-13/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_WT1{i,1} = bin_profile(raw_data,100);
end

% 5 min
for i = 1:21
    fid = fopen(strcat('/path/WT5_MG9-12-14/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_WT1{i,2} = bin_profile(raw_data,100);
end

% H1-TKO Pull of the replicates
all_chr_TKO1 = cell(21,2);

% 0 min
for i = 1:21
    fid = fopen(strcat('/path/TKO0_MG9-15-17/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_TKO1{i,1} = bin_profile(raw_data,100);
end

% 5 min
for i = 1:21
    fid = fopen(strcat('/path/TKO5_MG9-16-18/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_TKO1{i,2} = bin_profile(raw_data,100);
end

% Introduce information of the input gene list
% For each gene in the Input_genes.txt file
mygenes = tdfread('/path/Input_genes/Input_genes_20Kb.txt', "\t");
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
