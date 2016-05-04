% pVal = csvread('pvalue_adjacency_matrix.csv');
% [pL,pD] = build_Laplacian(pVal);

%take row 1 and make it into a list 
fileID = fopen('score_adjacency_matrix.csv', 'r');
 lineArr = fgetl(fileID);
 fclose(fileID);
 strArr = strrep(lineArr, ',','');
 genes = strsplit(strArr);
 %to look up specific index of list, do genes(i) where is index from 1 to max
 %genes(4400)


sVal = csvread('score_adjacency_matrix.csv', 1, 1);

M = abs(sVal);
n = length(sVal); d = 3; tmp = n^(-1/(d/2 + 2));
[L,~] = build_Laplacian(exp(M/tmp));

[U,D] = eig(L);
eigVals = diag(D);
[~,idx] = sort(eigVals,'descend');

k = 100;
cluster_idx = kmeans(U(:,idx(2:k)),k);
csvwrite(['cluster_idx',num2str(k),'.csv'],cluster_idx)









