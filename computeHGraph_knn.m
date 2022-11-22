function [L_hgraph,Dv,A] = computeHGraph_knn(X,knn,num_clu,flag)
% X: row denotes sample, column indicates features
% compute H with Knn
%H = Wtrim_knn(X,knn); 
% with louvain algorithm
H = Wtrim_louvain(X,num_clu,1);  %H = Wtrim_louvain(X,num_clu,2); 
W = Gaussion_kernel(X, H);
Dv = diag(sum(H*W,2));
temp = 1./(sum(H)+0.0001);
De = diag(temp);

if flag == 1
    % Zhou's method
    A = H*W*De*H';
elseif flag == 2
    % Rod
    A = H*W*H'-Dv;
    Dv = diag(sum(A));
end
L_hgraph = Dv - A;
% norm_LHGraph = Dv.^(-1/2)*L_hgraph*Dv.^(-1/2);

end


% Use Euclidean distance to compute the incidence matrix of hypergarph,
% including the node itself in KNN graph.
function H = Wtrim_knn(X,knn)
% X: row is sample
D = dist2(X, X); 
n_vertex = size(D,1);
[~,INDEX] = sort(D,2); % sort each row by ascend
H = zeros(n_vertex,n_vertex);
for i = 1:n_vertex
    H(i,INDEX(i, 2:knn+1)) = 1;
end

end


function H = Wtrim_louvain(X,num_clu,flag)
% X: row is sample
% Louvain to obtain initial cluteirng assignment
if flag == 1
    D = dist2(X,X);
    W = affinityMatrix(D,20); 
    A = Wtrim(W,50); 
    [clust,~,~] = getNCluster(A,num_clu,0,3,20); 
% DBSCAN 
elseif flag == 2
    epsilon = 0.5;
    MinPts = 1;
    clust = DBSCAN(X,epsilon,MinPts);
end

c = unique(clust);
num_clu = length(c);
n_vertex = size(X,1);
H = zeros(n_vertex, num_clu);
for i = 1:num_clu
    idx = clust == c(i);
    H(idx,i) = 1;
end

end

% construct hyperedge weight using heat kernel function
function W = Gaussion_kernel(X, H)
[~, n_HGedgs] = size(H);
W = zeros(n_HGedgs, n_HGedgs);
for i = 1:n_HGedgs
    idx = H(:,i)==1;
    X_new = X(idx,:);
    D = dist2(X_new, X_new); 
    delta = sum(sum(D))/(size(D,1)*(size(D,1)-1));
    y = exp(-D./(delta+1e-10));
    y0 = triu(y,1);
    W(i,i) = sum(sum(y0));
    
end
end



