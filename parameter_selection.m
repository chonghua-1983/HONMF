% parameter selection for OrthNMF algorithm
function [alpha, gamma, Inits] = parameter_selection(X1, X2, X3, A1, A2, A3, num_factor)
%num_clu = length(unique(label));
[~,H1] = nndsvd(A1,num_factor,0); [~,H2] = nndsvd(A2,num_factor,0);
[~,H3] = nndsvd(A3,num_factor,0);
G1 = eye(num_factor); G2 = eye(num_factor); G3 = eye(num_factor); 

Inits.H1 = H1'; Inits.H2 = H2'; Inits.H3 = H3';
%Inits.W1 = A1; Inits.W2 = A2; Inits.W3 = A3;
%Inits.D1 = diag(sum(A1)); Inits.D2 = diag(sum(A2)); Inits.D3 = diag(sum(A3)); 
knn = 12; flag = 1;
[L1,Dv,A] = computeHGraph_knn(X1',knn,num_factor,flag);
% Inits.D1 = Dv; Inits.W1 = A;
Inits.D1 = (L1+abs(L1))/2; Inits.W1 = (abs(L1)-L1)/2;
[L2,Dv,A] = computeHGraph_knn(X2',knn,num_factor,flag);
%Inits.D2 = Dv; Inits.W2 = A;
Inits.D2 = (L2+abs(L2))/2; Inits.W2 = (abs(L2)-L2)/2;
[L3,Dv,A] = computeHGraph_knn(X3',knn,num_factor,flag);
%Inits.D3 = Dv; Inits.W3 = A;
Inits.D3 = (L3+abs(L3))/2; Inits.W3 = (abs(L3)-L3)/2;
%L1 = Inits.D1 - Inits.W1; L2 = Inits.D2 - Inits.W2; L3 = Inits.D3 - Inits.W3;

Inits.G1 = G1; Inits.G2 = G2; Inits.G3 = G3;
Inits.S = SNF({A1,A2,A3},20);

err1 = norm(A1-H1'*G1*H1,'fro')^2; 
err2 = norm(A2-H2'*G2*H2,'fro')^2; 
err3 = norm(A3-H3'*G3*H3,'fro')^2; 

err4 = norm(Inits.S-H1'*H1,'fro')^2+norm(Inits.S-H2'*H2,'fro')^2+norm(Inits.S-H3'*H3,'fro')^2;
err5 = trace(H1*L1*H1'+H2*L2*H2'+H3*L3*H3');
alpha = (err1+err2+err3)/err4; 
gamma = abs((err1+err2+err3)/err5); 

% scaling parameter
%alpha = alpha*10;
%gamma = gamma/1000;
end
