function [W_F]=mergeW(W, sampleN, viewN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code is an implementation of the diffusion process in the paper:
%%% "Affinity Learning via a Diffusion Process for Subspace Clustering"
%%% By QILIN LI (li.qilin@postgrad.curtin.edu.au)
%%% Last Update 05/07/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WW = zeros(sampleN);
for i =1:viewN
    temp = W{i};
    WW = WW+temp;
end
W_F = WW;

end




