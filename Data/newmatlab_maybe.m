% pVal = csvread('am_test.csv');
% [pL,pD] = build_Laplacian(pVal);

%take row 1 and make it into a list 
fileID = fopen('am_test.csv', 'r');
 lineArr = fgetl(fileID);
 fclose(fileID);
 strArr = strrep(lineArr, ',','');
 genes = strsplit(strArr);
 %to look up specific index of list, do genes(i) where is index from 1 to max
 %genes(4400)


sVal = csvread('am_test.csv', 1, 1);

M = abs(sVal);
n = length(sVal); d = 3; tmp = n^(-1/(d/2 + 2));
[L,~] = build_Laplacian(exp(M/tmp));

[U,D] = eig(L);
eigVals = diag(D);
[~,idx] = sort(eigVals,'descend');

k = 3;
cluster_idx = kmeans(U(:,idx(2:k)),k);
csvwrite(['cluster_idx',num2str(k),'.csv'],cluster_idx);
     


genes;
size(genes,2)
size(cluster_idx,1)
cluster_idx;

output = cell(size(cluster_idx,1),2);
output;

size(cluster_idx,1);

cluster_idx(2);


for i = 1:size(cluster_idx,1)
   output(i,1) = genes(i+1);
   output(i,2) = cellstr(num2str(cluster_idx(i)));
   %output(i,2) = cluster_idx(i)
  
end
%sprintf(' %s %d', genes(i+1), cluster_idx(i))

output;

fileID = fopen('spec_output.txt', 'w');

formatSpec = '%s %s\r\n';

[nrows,ncols] = size(output);
for row = 1:nrows
    fprintf(fileID,formatSpec,output{row,:});
end

fclose(fileID);






