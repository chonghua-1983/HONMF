% select the neighborhoods of each node,(knn)
% K is the number of neighbors

function [W2] = Wtrim(W0,K)
[n,m] = size(W0);
% sort function, '2' is the row dimension,
% INDEX: the index value matrix of ranked W0
% W0: is the ordered matrix according to the INDEX
%[W0,INDEX]=sort(W0,2,'descend');
[W0,INDEX]=sort(W0,2,'descend');
W1 = zeros(n,n);
for i =1:m
%W1(i,INDEX(i,2:(K+1))) = W0(i,INDEX(i,2:(K+1)));
%W1(i,INDEX(i,1:K)) = W0(i,INDEX(i,1:K));
W1(i,INDEX(i,1:K)) = W0(i,1:K);
end
 W2 = (W1+W1')/2;