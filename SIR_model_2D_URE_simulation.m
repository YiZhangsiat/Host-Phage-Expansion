clear all;
close all;




A1_value=[2,3,4];
parfor i=1:length(A1_value)
    SIR_model_2D_URE_function(strcat(sprintf('initphage=1_GR=30_kap=2.5_PR=15_P0=0.05_Nk=5_Diff=%d.mat',A1_value(i))),A1_value(i));
end
