function [W_F]=mergeW(W, sampleN, viewN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code is an implementation of the mixture of diffusion affinity matrices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WW = zeros(sampleN);
for i =1:viewN
    temp = W{i};
    WW = WW+temp;
end
W_F = WW;

end




