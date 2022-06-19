% cd /media/cc/B/Josemi/TTseq_Feb2022/TTseq_scripts/Elongation_Rates_Fuchs/matlab
% Abre el archivo de escritura e imprime el header
outfid = fopen('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/Elongation_rate_5min_20220425_20Kb_size_Pull.txt', 'wt');
fprintf(outfid, "Gene_name\tWT_Pull\tTK0_Pull\n");
%fprintf(outfid, "Gene_name\tWT\tTKO\n");
%fprintf(outfid, "Gene_name\tWT\n");

% WT 1
all_chr_WT1 = cell(21,2);

for i = 1:21
    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/WT0_MG9-11-13/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_WT1{i,1} = bin_profile(raw_data,100);
end

for i = 1:21
    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/WT5_MG9-12-14/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_WT1{i,2} = bin_profile(raw_data,100);
end

% WT2
%all_chr_WT2 = cell(21,2);

%for i = 1:21
%    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate/1.2_Rate_calculation/MG9-13/profile_',num2str(i)));
%    raw_data = sparse(fscanf(fid, '%d'));
%    fclose(fid);
%    all_chr_WT2{i,1} = bin_profile(raw_data,100);
%end

%for i = 1:21
%    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate/1.2_Rate_calculation/MG9-14/profile_',num2str(i)));
%    raw_data = sparse(fscanf(fid, '%d'));
%   fclose(fid);
%    all_chr_WT2{i,2} = bin_profile(raw_data,100);
%end

% TKO1
all_chr_TKO1 = cell(21,2);

for i = 1:21
    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/TKO0_MG9-15-17/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_TKO1{i,1} = bin_profile(raw_data,100);
end

for i = 1:21
    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate_3/1.1_Rate_calculation/TKO5_MG9-16-18/profile_',num2str(i)));
    raw_data = sparse(fscanf(fid, '%d'));
    fclose(fid);
    all_chr_TKO1{i,2} = bin_profile(raw_data,100);
end

% TKO2
%all_chr_TKO2 = cell(21,2);

%for i = 1:21
%    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate/1.2_Rate_calculation/MG9-17/profile_',num2str(i)));
%    raw_data = sparse(fscanf(fid, '%d'));
%    fclose(fid);
%    all_chr_TKO2{i,1} = bin_profile(raw_data,100);
%end

%for i = 1:21
%    fid = fopen(strcat('/media/cc/B/Josemi/TTseq_Feb2022/TTseq_output/Elongation_rate/1.2_Rate_calculation/MG9-18/profile_',num2str(i)));
%    raw_data = sparse(fscanf(fid, '%d'));
%    fclose(fid);
%    all_chr_TKO2{i,2} = bin_profile(raw_data,100);
%end


% Para cada gen del archivo Input_genes.txt
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
    
    % El cuarto parámetro procede del archivo /media/cc/A/Josemi/PRUEBA_metlab/Elongation_rate/Alignment/RATIOS
    WT1bound=find_boundary_4sUDRBseq(gene_position, exons, all_chr_WT1, [1.17 1.35],100);
    %WT2bound=find_boundary_4sUDRBseq(gene_position, exons, all_chr_WT2, [1.00 1.62],100);
    TKO1bound=find_boundary_4sUDRBseq(gene_position, exons, all_chr_TKO1, [1.09 1.00],100);
    %TKO2bound=find_boundary_4sUDRBseq(gene_position, exons, all_chr_TKO2, [1.00 1.05],100);
    
    %fprintf(outfid, "%s\t%d\t%d\t%d\t%d\n", mygenes.Name(i,:), WT1bound, WT2bound, TKO1bound, TKO2bound);
    fprintf(outfid, "%s\t%d\t%d\n", mygenes.Name(i,:), WT1bound, TKO1bound);
    
    i = i+1;
end

% Cierra el archivo de escritura
fclose(outfid);
