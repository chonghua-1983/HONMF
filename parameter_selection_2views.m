% parameter selection for OrthNMF algorithm
function [alpha, gamma, Inits] = parameter_selection_2views(X1,X2,A1, A2, num_clu)
%num_clu = length(unique(label));
[~,H1] = nndsvd(A1,num_clu,0); [~,H2] = nndsvd(A2,num_clu,0);
G1 = eye(num_clu); G2 = eye(num_clu); 

Inits.H1 = H1'; Inits.H2 = H2'; 
% Inits.W1 = A1; Inits.W2 = A2; 
% Inits.D1 = diag(sum(A1)); Inits.D2 = diag(sum(A2)); 
% L1 = Inits.D1 - Inits.W1; L2 = Inits.D2 - Inits.W2; 

knn = 12; flag = 1;
[L1,Dv,A] = computeHGraph_knn(X1',knn,num_clu,flag);
Inits.D1 = (L1+abs(L1))/2; Inits.W1 = (abs(L1)-L1)/2;
[L2,Dv,A] = computeHGraph_knn(X2',knn,num_clu,flag);
Inits.D2 = (L2+abs(L2))/2; Inits.W2 = (abs(L2)-L2)/2;

Inits.G1 = G1; Inits.G2 = G2; 
Inits.S = SNF({A1,A2},20);

err1 = norm(A1-H1'*G1*H1,'fro')^2; 
err2 = norm(A2-H2'*G2*H2,'fro')^2; 

err4 = norm(Inits.S-H1'*H1,'fro')^2+norm(Inits.S-H2'*H2,'fro')^2;
err5 = trace(H1*L1*H1'+H2*L2*H2');
alpha = (err1+err2)/err4; 
gamma = abs((err1+err2)/err5); 

% scaling parameter
%alpha = alpha*10;
%gamma = gamma/10;
end