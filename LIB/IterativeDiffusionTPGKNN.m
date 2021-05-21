function [WW_G]=IterativeDiffusionTPGKNN(W,K, viewN)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code is an implementation of the diffusion process in the paper:
%%% "Affinity Learning via a Diffusion Process for Subspace Clustering"
%%% By QILIN LI (li.qilin@postgrad.curtin.edu.au)
%%% Last Update 05/07/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pre-processing of affinity matrix W
% d = sum(W, 2);
% D = diag(d + eps);
% W = W - diag(diag(W)) + D; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Normalization  %%%%%%%%%%%%%%%%%%
% W = W ./ repmat(sum(W, 2)+eps, 1, n);

for i = 1:viewN
    d = sum(W{i},2);
    D = diag(d + eps);
    W{i} = D^(-1/2)*W{i}*D^(-1/2); 
end

% d = sum(W,2);
% D = diag(d + eps);
% W = D^(-1/2)*W*D^(-1/2);      %%% Symmetric normalization is better
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = cell(viewN, 1);
for i =1:viewN
    S{i} = knnSparse(W{i}, K);
end
% S = knnSparse(W, K);          %%% Sparse is important for affinity matrix
WW_G = cell(viewN, 1);
for i=1:viewN
    switch i
        case 1
            WW = S{1};
            SS = (S{2}+S{3})/2;
        case 2
            WW = S{2};
            SS = (S{1}+S{3})/2;
        case 3
            WW = S{3};
            SS = (S{1}+S{2})/2;
    end
    maxIter = 50;
    epsilon = 1e-2;        %%% convergence threshold
    for t = 1:maxIter
        temp = SS*WW*SS' + eye(length(WW));
        if norm(temp-WW,'fro') < epsilon, break; end  
        WW = temp;   
    end
    WW_G{i} = knnSparse(WW, K); 
end
% WW = S;    
% maxIter = 50;
% epsilon = 1e-2;        %%% convergence threshold
% for t = 1:maxIter
%     temp = S*WW*S' + eye(length(WW));
%     if norm(temp-WW,'fro') < epsilon, break; end  
%     WW = temp;   
% end
% WW = knnSparse(WW, K);