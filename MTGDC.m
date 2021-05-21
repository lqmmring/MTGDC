function [W,WW_Mer,out,RES,Con]=MTGDC(X,pr,label)
ViewN = 3;
kmeansK = length(unique(label));
TotalSampleNo=length(label);
TempvData=X;
NorTempvData=NormalizeFea(double(TempvData));
[tempN,~] = size(TempvData);
if(pr==1)
    k=20;
else
    k=ceil(tempN*pr);
end

tic
disp('Affinity learning......');
if(pr==1)
   [W,~] = adaptiveGaussian(NorTempvData, k, ViewN);
   WW=IterativeDiffusionTPGKNN(W,k, ViewN);
   WW_Mer=mergeW(WW,TotalSampleNo,ViewN);
else
   [~,W] = adaptiveGaussian(NorTempvData, k, ViewN);
   WW=IterativeDiffusionTPGKNN(W,k, ViewN);
   WW_Mer=mergeW(WW,TotalSampleNo,ViewN);
end
toc

tic
disp('Spectral clustering......');
out = SpectralClustering(WW_Mer,kmeansK);
toc

% [8: ACC MIhat Purity ARI F-score Precision Recall Contingency];
[result,Con] = ClusteringMeasure(label, out');  
nmi = Cal_NMI(label, out');
RES=[result, nmi];
end