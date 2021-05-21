function [W,W_knn]=adaptiveGaussian(data,K,viewN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code is an implementation of adaptive Gaussian affinity used in:
%%% "Affinity Learning via a Diffusion Process for Subspace Clustering"
%%% Note that the self-affinity is defined, and it's not simply 0 or 1
%%% By Qiaoming Liu (cslqm@hit.edu.cn)
%%% Last Update 07/12/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = size(data, 1);
W =cell(viewN, 1);
W_knn =cell(viewN, 1);

for i=1:viewN
    switch i
        case 1
            D = EuDist2(data);
            
        case 2
            distance = pdist(data,'correlation');
            D = squareform(distance);
            
        case 3
            distance = pdist(data,'spearman');
            D = squareform(distance);
    end
    
    D = D - diag(diag(D));   %%% Zero distance to itself
    [T, idx] = sort(D, 2);
    
    temp_W = zeros(n,n);
    for q = 1:n
        for j = 1:n
            sigma = mean([T(q,2:K+1), T(j,2:K+1)]);
            temp_W(q,j) = normpdf(D(q,j), 0, 0.35*sigma);
        end
    end
    W{i} = ( temp_W+ temp_W')/2;
    P = zeros(n, n);
    for z = 1:n
        P(z, idx(z,1:K+1)) = 1;
    end

    temp_W_knn = temp_W.*P;
    W_knn{i} = (temp_W_knn+temp_W_knn')/2;
end