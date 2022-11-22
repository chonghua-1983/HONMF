function [H1,H2,H3,S,G1,G2,G3,objs,iter] = OrthNMF(A1,A2,A3,alpha, gamma,Inits)
% Orthogonal NMF for multi-view microbiome data
% A1: sample-sample similarity matrix bacterial
% A2, A3: fungi, vircus
% H1,H2,H3£ºlow-dimensional sample representation matrices
% S: consensus sample similarity matrix
% G1: digonal matrix
H1 = Inits.H1; H2 = Inits.H2; H3 = Inits.H3; 
G1 = Inits.G1; G2 = Inits.G2; G3 = Inits.G3;
D1 = Inits.D1; D2 = Inits.D2; D3 = Inits.D3;  
W1 = Inits.W1; W2 = Inits.W2; W3 = Inits.W3; S = Inits.S;
L1 = D1-W1; L2 = D2-W2; L3 = D3-W3;
n = size(S,1);

obj_old = 1; stop_rule = 2; 
clear Inits

disp('iteration starts!');
Maxiter = 2000; objs = zeros(Maxiter,1); yita = 10; %yita = 10
for iter = 1:Maxiter
    % update H1 H2 H3
    H1tH1 = H1'*H1; H2tH2 = H2'*H2; H3tH3 = H3'*H3; 
    H1 = H1.*(A1*H1*G1 + yita*H1 + alpha*S'*H1 + 0.5*gamma*W1*H1)./(H1*G1*H1tH1*G1 + 0.5*gamma*D1*H1 + (alpha+yita)*H1*H1tH1);
    H2 = H2.*(A2*H2*G2 + yita*H2 + alpha*S'*H2 + 0.5*gamma*W2*H2)./(H2*G2*H2tH2*G2 + 0.5*gamma*D2*H2 + (alpha+yita)*H2*H2tH2);
    H3 = H3.*(A3*H3*G3 + yita*H3 + alpha*S'*H3 + 0.5*gamma*W3*H3)./(H3*G3*H3tH3*G3 + 0.5*gamma*D3*H3 + (alpha+yita)*H3*H3tH3);
    % update G1 G2 G3
    G1 = G1.*(H1'*A1*H1./((H1'*H1)*G1*(H1'*H1)));
    G2 = G2.*(H2'*A2*H2./((H2'*H2)*G2*(H2'*H2)));
    G3 = G3.*(H3'*A3*H3./((H3'*H3)*G3*(H3'*H3)));
    % update S
    H1H1t = H1*H1'; H2H2t = H2*H2'; H3H3t = H3*H3'; 
    tmp = alpha*(H1H1t + H2H2t + H3H3t);
    sum_S = repmat(sum(S,2),1,n); 
    S = S.*(tmp + ones(n))./(3*alpha*S + sum_S);
    
    if stop_rule == 2
        obj = compute_obj(A1,A2,A3,S,H1,H2,H3,G1,G2,G3,L1,L2,L3,alpha,gamma,yita);
        objs(iter,1) = obj;
        if (abs(obj_old-obj)/obj_old < 10^(-5) && iter > 1) || iter == Maxiter
            disp('converged!');
            break;
        end
        obj_old = obj;
    end
    if mod(iter,50) == 0
        disp(['number of iteration:',num2str(iter),'  obj:',num2str(obj)]);
    end
    
end

S = (S+S')/2;
end

function obj = compute_obj(A1,A2,A3,S,H1,H2,H3,G1,G2,G3,L1,L2,L3,alpha,gamma,yita)
  n = size(S,1); k = size(H1,2);
  error_1 = norm(A1-H1*G1*H1','fro')^2+norm(A2-H2*G2*H2','fro')^2+norm(A3-H3*G3*H3','fro')^2;
  error_2 = alpha*(norm(S-H1*H1','fro')^2+norm(S-H2*H2','fro')^2+norm(S-H3*H3','fro')^2);
  error_3 = gamma*(trace(H1'*L1*H1+H2'*L2*H2+H3'*L3*H3));
  error_4 = yita*(norm(H1'*H1-eye(k),'fro')^2 + norm(H2'*H2-eye(k),'fro')^2 +norm(H3'*H3-eye(k),'fro')^2);
  error_5 = (norm(sum(S,2) - ones(n,1),'fro')^2);
  obj = error_1 + error_2 + error_3+ error_4 + error_5;
end

