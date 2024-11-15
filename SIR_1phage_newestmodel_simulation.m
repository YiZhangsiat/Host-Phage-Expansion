clear all;
close all;

%  SIR_1phage_Nutr_Phag_HillFunct_function('NutrPhag_Hill_init=1_A1=51.mat',51);




A1_value=[96,97,98,99];
parfor i=1:length(A1_value)
    SIR_1phage_newestmodel_function(strcat(sprintf('newestmodel_constInfection_initphage=1E-1_Dp=0_attr=100_chi=300_phagePR=%d.mat',A1_value(i))),A1_value(i));
end


% % % A1_value =[28 30:1:50];
% % A1_value =[ 35 ];
% % 
% % % A2_value = [28 30:1:50];
% % A2_value =[ 35.1 ];
% % 
% 
% for i = 1:length(A1_value)
%     
%  for j = 1:length(A2_value)
%    SIR_1phage_Nutr_Phag_HillFunct_function(strcat(sprintf('NutrPhag_Hill_init=1_A1=%d',A1_value(i)),sprintf('_A2=%.2e.mat',A2_value(j))),A1_value(i),A2_value(j))
% 
%  end
%  
% end
