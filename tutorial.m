% This is a demo for implementing HONMF on multi-omics microbiome dataset
load('data1.mat')
load('data1_label.mat')
addpath('external/')
num_clu = length(unique(label));
n = size(X1,2);

D1 = dist2(X1',X1'); A1 = affinityMatrix(D1,20); 
D2 = dist2(X2',X2'); A2 = affinityMatrix(D2,20); 
D3 = dist2(X3',X3'); A3 = affinityMatrix(D3,20); 
[alpha, gamma, Inits] = parameter_selection(X1, X2, X3, A1, A2, A3, num_clu);

[H1,H2,H3,S,G1,G2,G3,objs,iter] = OrthNMF(A1,A2,A3,alpha, gamma,Inits);
A = Wtrim(S,29);  % half of the number of samples
[clust,~,~] = getNCluster(A,num_clu,0,3,20); 
if length(unique(clust))== num_clu
   [ac,ari, nmi_value, ~] = CalcMetrics(label, clust); 
else
   [~, clust, ~] = SpectralClustering(A, num_clu); 
   [ac,ari, nmi_value, ~] = CalcMetrics(label, clust); 
end
% save clusteirng results
csvwrite("data1_A.csv",A)
csvwrite("data1_S.csv",S)
csvwrite("data1_clust.csv",clust)

% UMAP visualization
addpath('umapFileExchange/umap') 
addpath('umapFileExchange/util')

term = '.'; 
label_name = load('data1.mat','label');
label_name = label_name.label;
%label_name = strrep(label_name,'_','\_');
colors = generateColors(max(length(unique(label_name))),length(unique(term)));

D_orthnnmf = 1-S; D_orthnnmf = D_orthnnmf-diag(diag(D_orthnnmf));
[reduction_jsnmf, ~, ~, ~] = run_umap(D_orthnnmf,'metric','precomputed'); %,'min_dist',0.68,'n_neighbors',12
% reduction_jsnmf = tsne_p(S);
gscatter(reduction_jsnmf(:,1),reduction_jsnmf(:,2),label_name,colors,[],10); 
set(gca,'xtick',[],'ytick',[]);
title('UMAP for JSNMF')

legend('Location','westoutside','Box','off','FontSize',9.5); 
legendmarkeradjust(16)
legend('boxoff') 

% display clustering
displayClusters(S,clust,1); % flag = 1

% parameter sensivity analysis
alpha_set = [alpha/10,alpha/5,alpha/2,alpha, alpha*2, alpha*5, alpha*10];
c = 0; len = length(alpha_set);
results_alpha = zeros(len, 3);
for i = 1:length(alpha_set)
    c = c +1;
    results_alpha(c, 1) = alpha_set(i);
    [H1,H2,H3,S,G1,G2,G3,objs,iter] = OrthNMF(A1,A2,A3,alpha_set(i), gamma,Inits);
    A = Wtrim(S,29); 
    [clust,~,~] = getNCluster(A,num_clu,0,3,20); 
    if length(unique(clust))== num_clu
       [ac,ari, nmi_value, ~] = CalcMetrics(label, clust); 
    else
       [~, clust, ~] = SpectralClustering(A, num_clu); 
       [ac,ari, nmi_value, ~] = CalcMetrics(label, clust); 
    end
    results_alpha(c, 2) = nmi_value;
    results_alpha(c, 3) = ari;
end

gamma_set = [gamma/10, gamma/5,gamma/2, gamma, gamma*2, gamma*5,gamma*10];
c = 0; len = length(gamma_set);
results_gamma = zeros(len, 3);
for i = 1:length(gamma_set)
    c = c +1;
    results_gamma(c, 1) = gamma_set(i);
    [H1,H2,H3,S,G1,G2,G3,objs,iter] = OrthNMF(A1,A2,A3,alpha, gamma_set(i),Inits);
    A = Wtrim(S,29); 
    [clust,~,~] = getNCluster(A,num_clu,0,3,20); 
    if length(unique(clust))== num_clu
       [ac,ari, nmi_value, ~] = CalcMetrics(label, clust); 
    else
       [~, clust, ~] = SpectralClustering(A, num_clu); 
       [ac,ari, nmi_value, ~] = CalcMetrics(label, clust); 
    end
    results_gamma(c, 2) = nmi_value;
    results_gamma(c, 3) = ari;
end

% Plot the curve performance of HONMF vs alpha
log_alpha = alpha_set;
nmi2= results_alpha(:,2); ari2 = results_alpha(:,3);
plot(1:length(log_alpha),nmi2,'-*b',1:length(log_alpha),ari2,'-ok','MarkerSize',4)
axis([0.75 7.25 0.2 0.9])

lgd = legend('NMI', 'ARI');
set(lgd,'FontName','Times New Roman','FontSize',9,'FontWeight','normal','Location','NorthEast')
%axis([0.8 7.2 0.38 0.87]); % alter the range of X and Y axis

set(gca,'Xtick',1:length(log_alpha),'Ytick',[0.3 0.5 0.7 0.9], 'fontsize',9,'fontname','Times New Roman','FontWeight','bold');
set(gca,'Xticklabel',{'\alpha*/10','\alpha*/5','\alpha*/2','\alpha*','2\alpha*','5\alpha*','10\alpha*'}, 'Yticklabel',{'0.3','0.5','0.7','0.9'},'FontWeight','bold');

%xlabel('log \alpha','FontSize',10,'FontWeight','bold');
ylabel('Clustering performnce','FontSize',9,'FontWeight','bold');
title('The performance of HONMF w.r.t \alpha');

% Plot the performance curve of HONMF vs gamma
log_gamma = gamma_set;
nmi2= results_gamma(:,2); ari2 = results_gamma(:,3);
plot(1:length(log_gamma),nmi2,'-*b',1:length(log_gamma),ari2,'-ok','MarkerSize',4)
axis([0.75 7.25 0.2 0.9])

lgd = legend('NMI', 'ARI');
set(lgd,'FontName','Times New Roman','FontSize',9,'FontWeight','normal','Location','NorthEast')
%axis([0.8 7.2 0.38 0.87]); % alter the range of X and Y axis

set(gca,'Xtick',1:length(log_gamma),'Ytick',[0.3 0.5 0.7 0.9], 'fontsize',9,'fontname','Times New Roman','FontWeight','bold');
set(gca,'Xticklabel',{'\gamma*/10','\gamma*/5','\gamma*/2','\gamma*','2\gamma*','5\gamma*','10\gamma*'}, 'Yticklabel',{'0.3','0.5','0.7','0.9'},'FontWeight','bold');

%xlabel('log \alpha','FontSize',10,'FontWeight','bold');
ylabel('Clustering performnce','FontSize',9,'FontWeight','bold');
title('The performance of HONMF w.r.t \gamma');

% Runing HONMF on 2 views microbiome multi-omics data
[alpha, gamma, Inits] = parameter_selection_2views(X1', X2', A1, A2, num_clu);
[H1,H2,S,G1,G2,objs,iter] = OrthNMF_2views(A1,A2,alpha, gamma,Inits);

term = '.'; 
colors = generateColors(max(num_clu,length(unique(term))));

D_orthnnmf = 1-S; D_orthnnmf = D_orthnnmf-diag(diag(D_orthnnmf)); D_orthnnmf(D_orthnnmf<0) = 0;
[reduction_jsnmf, ~, ~, ~] = run_umap(D_orthnnmf,'metric','precomputed','min_dist',0.68,'n_neighbors',12); %,'min_dist',0.68,'n_neighbors',12
gscatter(reduction_jsnmf(:,1),reduction_jsnmf(:,2),clust,colors,[],4); 
set(gca,'xtick',[],'ytick',[]);
title('UMAP for JSNMF')

legend('Location','westoutside','Box','off','FontSize',9.5); 
legendmarkeradjust(16)
legend('boxoff') 
