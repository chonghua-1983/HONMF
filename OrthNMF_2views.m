function [H1,H2,S,G1,G2,objs,iter] = OrthNMF_2views(A1,A2,alpha, gamma,Inits)
% Orthogonal NMF for multi-view microbiome data
% A1: sample-sample similarity matrix bacterial
% A2, A3: fungi, vircus
% H1,H2,H3£ºlow-dimensional sample representation matrices
% S: consensus sample similarity matrix
% G1: digonal matrix
H1 = Inits.H1; H2 = Inits.H2; 
G1 = Inits.G1; G2 = Inits.G2;
D1 = Inits.D1; D2 = Inits.D2; 
W1 = Inits.W1; W2 = Inits.W2; S = Inits.S;
L1 = D1-W1; L2 = D2-W2; 
n = size(S,1);

obj_old = 1; stop_rule = 2; 
clear Inits

disp('iteration starts!');
Maxiter = 500; objs = zeros(Maxiter,1); yita = 2; %yita = 2
for iter = 1:Maxiter
    % update H1 H2
    H1tH1 = H1'*H1; H2tH2 = H2'*H2; 
    H1 = H1.*(A1*H1*G1 + H1 + alpha*S'*H1 + 0.5*gamma*W1*H1)./(H1*G1*H1tH1*G1 + 0.5*gamma*D1*H1 + (alpha+yita)*H1*H1tH1);
    H2 = H2.*(A2*H2*G2 + H2 + alpha*S'*H2 + 0.5*gamma*W2*H2)./(H2*G2*H2tH2*G2 + 0.5*gamma*D2*H2 + (alpha+yita)*H2*H2tH2);
    % update G1 G2
    G1 = G1.*(H1'*A1*H1./((H1'*H1)*G1*(H1'*H1)));
    G2 = G2.*(H2'*A2*H2./((H2'*H2)*G2*(H2'*H2)));
    % update S
    H1H1t = H1*H1'; H2H2t = H2*H2';
    tmp = alpha*(H1H1t + H2H2t);
    sum_S = repmat(sum(S,2),1,n); 
    S = S.*(tmp + ones(n))./(2*alpha*S + sum_S);
    
    if stop_rule == 2
        obj = compute_obj(A1,A2,S,H1,H2,G1,G2,L1,L2,alpha,gamma,yita);
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

function obj = compute_obj(A1,A2,S,H1,H2,G1,G2,L1,L2,alpha,gamma,yita)
  n = size(S,1); k = size(H1,2);
  error_1 = norm(A1-H1*G1*H1','fro')^2+norm(A2-H2*G2*H2','fro')^2;
  error_2 = alpha*(norm(S-H1*H1','fro')^2+norm(S-H2*H2','fro')^2);
  error_3 = gamma*(trace(H1'*L1*H1+H2'*L2*H2));
  error_4 = yita*(norm(H1'*H1-eye(k),'fro')^2 + norm(H2'*H2-eye(k),'fro')^2);
  error_5 = (norm(sum(S,2) - ones(n,1),'fro')^2);
  obj = error_1 + error_2 + error_3+ error_4 + error_5;
end

