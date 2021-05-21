clc;
clear;
close all;
addpath('MeasureTools');
addpath('LIB');
% compile_func(0);
load('Data/Data_Buettner.mat'); % load data
X=in_X;
pr = 1;
label=true_labs;
[W,WW_Mer,out,RES,Con]=MTGDC(X,pr,label);