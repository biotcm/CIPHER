%% Cipher
% Written by: Xuebing Wu, 2008
% Revised by: Aidi Tan, 2013

%% Initialization
clc
clear

inner_pheno_pheno_similarity = load('inner_phenotype_similarity.mat');
fname = fieldnames(inner_pheno_pheno_similarity);
inner_pheno_pheno_similarity = inner_pheno_pheno_similarity.(fname{1});
phenoNum = length(inner_pheno_pheno_similarity);

clear fname
%% Build PPI Network
disp('Computing The Shortest Distance Between Protein Nodes...');
tic

% Load
ppi_matrix_data = load('inner_ppi.txt');
geneNum = max(max(ppi_matrix_data));
ppi_matrix = sparse(ppi_matrix_data(:,1), ppi_matrix_data(:,2), 1, geneNum, geneNum);
ppi_matrix = ppi_matrix + ppi_matrix';

% Compute PPI_Shortest_Dist
Protein_Shortest_Distance = graphallshortestpaths(ppi_matrix, 'Directed', 'false');

toc
clear ppi_matrix_data
%% Calculate Phenotype Gene Closeness
disp('Calculate Phenotype Gene Closeness...');
tic
pheno_gene_relation = cell(1, phenoNum);

% Load
fid = fopen('inner_phenotype_gene_relation.txt'); %Get file id
line = fgetl(fid);
while ischar(line) && ~isempty(line)
    chArray = regexp(line,'\t','split'); %Get target array for this drug
    phenoIndex = str2double(chArray(1));
    gNum = size(chArray,2) - 1;
    Array = zeros(1, gNum);
    for i = 1:gNum
        Array(i) = str2double(chArray{i+1});
    end
    pheno_gene_relation{phenoIndex} = Array;

    % Next line
    line = fgetl(fid);
end
fclose(fid);
clear fid line chArray Array GeneNum ans i

% Compute Gene2Phenotype Closeness
gene2phenotype_closeness = zeros(geneNum, phenoNum);
for phenoIndex = 1:phenoNum
    Array = pheno_gene_relation{phenoIndex};
    if isempty(Array)
        gene2phenotype_closeness(:, phenoIndex) = 0;
    else
        for geneIndex = 1:geneNum
            gene2phenotype_closeness(geneIndex, phenoIndex) = sum(exp(-(Protein_Shortest_Distance(geneIndex,Array)).^2));
        end
    end
end

toc
clear geneIndex phenoIndex ans array gNum
%% Compute Phenotype Gene Score
disp('Compute Phenotype Gene Score...')
tic
Phenotype_Gene_Score_Matrix = corr(inner_pheno_pheno_similarity, gene2phenotype_closeness');
toc
nan_row = [];
for phenoIndex = 1:phenoNum
    if sum(isnan(Phenotype_Gene_Score_Matrix(phenoIndex,:))) == geneNum
        nan_row = [nan_row, phenoIndex];
    end
end
Phenotype_Gene_Score_Matrix(isnan(Phenotype_Gene_Score_Matrix)) = -1;
Phenotype_Gene_Score_Matrix(nan_row, :) = nan;
save('inner_pheno_gene_score.mat','Phenotype_Gene_Score_Matrix');
dlmwrite('inner_pheno_gene_score.txt', Phenotype_Gene_Score_Matrix,'\t');
