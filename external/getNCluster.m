function [clust,Q,i] = getNCluster(W,n_cluster,range_min,range_max,max_iter)
% Louvain algorithm with resolution parameter
% input:
% W: similarity matrix
% n_cluster: number of 
% range_min = 0; range_max = 3;
% max_iter: number of interates

for i =1:max_iter
    solution = range_min+(range_max-range_min)/2;
    [M,Q]=community_louvain(W,solution);
    this_cluster=length(unique(M));
    if this_cluster > n_cluster
        range_max = solution;
    elseif this_cluster < n_cluster
        range_min = solution;
    else
        clust = M;
        break;
    end
    if i == max_iter
        clust = M;
        break;
    end
end

end




