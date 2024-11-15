
clear all;
close all;


TimeT=15:32;
PR =[85];


for m=1:length(TimeT)
 
 filename=strcat('newestmodel_constInfection_initphage=1E-1_Dp=0_attr=100_chi=300_phagePR=',int2str(PR),'_',int2str(TimeT(m)));
 filename1=strcat('newestmodel_constInfection_initphage=1E-1_Dp=0_attr=100_chi=300_phagePR=',int2str(PR),'_',int2str(TimeT(m)));
load(strcat(filename1,'.mat'));

     %%  
     Total_cell=Cell_den_S+Cell_den_I+Cell_den_R;
     infected_cell=Cell_den_I+Cell_den_R;
     

     [Nx,Ny]=size(Total_cell);

     Nx0=round(Nx/2); 
     Ny0=round(Ny/2);
     

     %% 
     
         [~,FW_exp_loca]=findpeaks(Total_cell(:,Ny0));

        FW_exp_l=FW_exp_loca(end);
      cell_infected_peak_location(m)=FW_exp_l*dx/10000; %cm
      cell_infected_peak_value(m)=Total_cell(FW_exp_l,Ny0);
     
     
   %% 

       [~,FW_con_loca]=findpeaks(Total_cell(Nx0,:));  
       FW_con_l=FW_con_loca(end);

     cell_con_peak_location(m)=FW_con_l*dx/10000; %cm
     cell_con_peak_value(m)=Total_cell(Nx0,FW_con_l);

      
end


  %% 
  real_time_h=TimeT;
cell_infected_peak_time=real_time_h;
cell_infected_peak_time(isnan(cell_infected_peak_location))=[];
cell_infected_peak_location(isnan(cell_infected_peak_location))=[];
MR_cell_infected_para=polyfit(cell_infected_peak_time,cell_infected_peak_location,1);
MR_cell_infected=MR_cell_infected_para(1)*cell_infected_peak_time+MR_cell_infected_para(2);
cell_infected_MR_V=MR_cell_infected_para(1);
% R2_Cell_infected_line(i,:)=norm(MR_cell_infected-mean(cell_infected_peak_location))^2/norm(cell_infected_peak_location-mean(cell_infected_peak_location))^2;
%% 
cell_con_peak_time=real_time_h;
cell_con_peak_time(isnan(cell_con_peak_location))=[];
cell_con_peak_location(isnan(cell_con_peak_location))=[];
MR_cell_con_para=polyfit(cell_con_peak_time,cell_con_peak_location,1);
MR_cell_con=MR_cell_con_para(1)*cell_con_peak_time+MR_cell_con_para(2);
cell_con_MR_V=MR_cell_con_para(1);
% R2_Cell_con_line(i,:)=norm(MR_cell_con-mean(cell_con_peak_location))^2/norm(cell_con_peak_location-mean(cell_con_peak_location))^2;


%%  
Cell_Exp_peak_infor=figure('position',[100 100 500 500]);
subplot(2,1,1);plot(real_time_h,cell_infected_peak_value); 
text(real_time_h(5),cell_infected_peak_value(5),sprintf('peak value=%.3f',max(cell_infected_peak_value)),'FontSize',30);

xlabel('time');
ylabel('Cell infected peak value');
title('cell infected front peak')
subplot(2,1,2);plot(cell_infected_peak_time,cell_infected_peak_location); 
hold on;plot(cell_infected_peak_time,MR_cell_infected,'r');
text(real_time_h(5),cell_infected_peak_location(5),sprintf('MR=%.3f',cell_infected_MR_V),'FontSize',30);

xlabel('time');
ylabel('Cell infected peak location');

saveas(Cell_Exp_peak_infor,strcat(filename,'_Exp_MR'),'fig');
saveas(Cell_Exp_peak_infor,strcat(filename,'_Exp_MR'),'png');

%% 
Cell_con_peak_infor=figure('position',[100 100 500 500]);
subplot(2,1,1);plot(real_time_h,cell_con_peak_value); 
text(real_time_h(5),cell_con_peak_value(5),sprintf('peak value=%.3f',max(cell_con_peak_value)),'FontSize',30);

xlabel('time');
ylabel('cell con peak_value');
title('con cell front peak')
subplot(2,1,2);plot(cell_con_peak_time,cell_con_peak_location); 
hold on;plot(cell_con_peak_time,MR_cell_con,'r');
text(real_time_h(5),cell_con_peak_location(end),sprintf('MR=%.3f',cell_con_MR_V),'FontSize',30);

xlabel('time');
ylabel('Cell con peak location');
saveas(Cell_con_peak_infor,strcat(filename,'_con_MR'),'fig');
saveas(Cell_con_peak_infor,strcat(filename,'_con_MR'),'png');

