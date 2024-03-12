clear all;
close all;



% 
% SIR_model_2D_simplifyV4_function('L=9cm_chi=370_alpha=50.mat',43);
A1_value=[7,9,11,13,15];
parfor i=1:length(A1_value)
    SIR_model_2D_simplifyV4_function(strcat(sprintf('C_dx=100_R=0.8_tar35Simu_initP=1E-3_diff=30_K=3.5_PR=10_chi=180_eta=%d.mat',A1_value(i))),A1_value(i));
end
