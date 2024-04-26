clear all;
close all;


TimeT=220:250;
Diff_value =[2];%[60 90 120 150 180];

for i=1:length(Diff_value)
for m=1:length(TimeT)

    filename=strcat('initphage=1_GR=30_kap=2p5_PR=15_P0=0p05_Nk=5_Diff=',int2str(Diff_value(i)),'_',int2str(TimeT(m)));
    filename1=strcat('initphage=1_GR=30_kap=2.5_PR=15_P0=0.05_Nk=5_Diff=',int2str(Diff_value(i)),'_',int2str(TimeT(m)));

     load(strcat(filename1,'.mat'));
     %%  
     Total_cell=Cell_den_S+Cell_den_I+Cell_den_R;
     infected_cell=Cell_den_I+Cell_den_R;
     
     ref_value=20;

  
     
     [~,FW_loca]=min(abs(Total_cell(ref_value:end,50)-1));  %计算fisherwave上total density为0.9的位置
      FW_l=FW_loca(end)+ref_value-1;
          
     %
     infect_ref_value=0.0001;
%      LY_C=round(round(LY/dx)/2);
%          [~,FW_Cell_Infect_half_loc_B]=min(abs(infected_cell(FW_l,1:LY_C)- infect_ref_value));
%          [~,FW_Cell_Infect_half_loc_E]=min(abs(infected_cell(FW_l,LY_C:end)- infect_ref_value));
%          
%          FW_location_x(m)=FW_l;
%          FW_location_y_left(m)=FW_Cell_Infect_half_loc_B;
%          FW_location_y_right(m)=FW_Cell_Infect_half_loc_E+LY_C-1;
% 
%          FW_Infection_range(m)=(FW_location_y_right(m)-FW_location_y_left(m))*dx/10000;      
  
     
     
     [FWpeak_v,FWpeak_l]=findpeaks(infected_cell(FW_l,:));
     
     
      if ~isnan(length(FWpeak_l))
          
         [~,FW_Cell_Infect_half_loc_B]=min(abs(infected_cell(FW_l,1:FWpeak_l(1))-infect_ref_value));
         [~,FW_Cell_Infect_half_loc_E]=min(abs(infected_cell(FW_l,FWpeak_l(end):end)-infect_ref_value));
         
         FW_location_x(m)=FW_l;
         FW_location_y_left(m)=FW_Cell_Infect_half_loc_B;
         FW_location_y_right(m)=FW_Cell_Infect_half_loc_E+FWpeak_l(end)-1;

         FW_Infection_range(m)=(FW_location_y_right(m)-FW_location_y_left(m))*dx/10000; 
      
      end

     %%  Exp_MR 
         fixed_exp_V=round(LY/dy)/2-50;
     
         [~,FW_exp_loca]=min(abs(Total_cell(ref_value:end,fixed_exp_V)-0.5));  %计算fisherwave上total density为0.9的位置
         FW_exp_l=FW_exp_loca(end)+ref_value-1;

         cell_infected_peak_location(m)=FW_exp_l*dx/10000; %cm
         cell_infected_peak_value(m)=Total_cell(FW_exp_l,fixed_exp_V);
     
     
   %% Con_MR
     FW_con_l=FW_l;
     cell_con_peak_location(m)=FW_con_l*dx/10000; %cm
     cell_con_peak_value(m)=Total_cell(FW_con_l,100);
      
end
end

  %% migration velocity Cell 
  real_time_h=TimeT;
cell_infected_peak_time=real_time_h;
cell_infected_peak_time(isnan(cell_infected_peak_location))=[];
cell_infected_peak_location(isnan(cell_infected_peak_location))=[];
MR_cell_infected_para=polyfit(cell_infected_peak_time,cell_infected_peak_location,1);
MR_cell_infected=MR_cell_infected_para(1)*cell_infected_peak_time+MR_cell_infected_para(2);
cell_infected_MR_V=MR_cell_infected_para(1);
% R2_Cell_infected_line(i,:)=norm(MR_cell_infected-mean(cell_infected_peak_location))^2/norm(cell_infected_peak_location-mean(cell_infected_peak_location))^2;
%% migration velocity Cell 
cell_con_peak_time=real_time_h;
cell_con_peak_time(isnan(cell_con_peak_location))=[];
cell_con_peak_location(isnan(cell_con_peak_location))=[];
MR_cell_con_para=polyfit(cell_con_peak_time,cell_con_peak_location,1);
MR_cell_con=MR_cell_con_para(1)*cell_con_peak_time+MR_cell_con_para(2);
cell_con_MR_V=MR_cell_con_para(1);
% R2_Cell_con_line(i,:)=norm(MR_cell_con-mean(cell_con_peak_location))^2/norm(cell_con_peak_location-mean(cell_con_peak_location))^2;

%% infection speed  V//
cell_infection_time=real_time_h;
Infection_speed_para=polyfit(cell_infection_time,FW_Infection_range,1);
Infection_speed_line=Infection_speed_para(1)*cell_infection_time+Infection_speed_para(2);
Infection_speed_V=Infection_speed_para(1);



%%  
Cell_Exp_peak_infor=figure('position',[100 100 500 500]);
subplot(2,1,1);plot(real_time_h,cell_infected_peak_value); 
text(real_time_h(end-10),cell_infected_peak_value(end-10),sprintf('peak value=%.3f',max(cell_infected_peak_value)),'FontSize',30);
% legend('Cell 1 peak value');
xlabel('time');
ylabel('Cell infected peak value');
title('cell infected front peak')
subplot(2,1,2);plot(cell_infected_peak_time,cell_infected_peak_location); 
hold on;plot(cell_infected_peak_time,MR_cell_infected,'r');
text(real_time_h(end-5),cell_infected_peak_location(end-5),sprintf('MR=%.3f',cell_infected_MR_V),'FontSize',30);
% legend('Cell 1 major peak location','Cell 1 line');
xlabel('time');
ylabel('Cell infected peak location');

saveas(Cell_Exp_peak_infor,strcat(filename,'_Exp_MR'),'fig');
saveas(Cell_Exp_peak_infor,strcat(filename,'_Exp_MR'),'png');

%% 
Cell_con_peak_infor=figure('position',[100 100 500 500]);
subplot(2,1,1);plot(real_time_h,cell_con_peak_value); 
text(real_time_h(end-10),cell_con_peak_value(end-10),sprintf('peak value=%.3f',max(cell_con_peak_value)),'FontSize',30);
% legend('Cell 1 peak value');
xlabel('time');
ylabel('cell con peak_value');
title('con cell front peak')
subplot(2,1,2);plot(cell_con_peak_time,cell_con_peak_location); 
hold on;plot(cell_con_peak_time,MR_cell_con,'r');
text(real_time_h(end-10),cell_con_peak_location(end),sprintf('MR=%.3f',cell_con_MR_V),'FontSize',30);
% legend('Cell 1 major peak location','Cell 1 line');
xlabel('time');
ylabel('Cell con peak location');
saveas(Cell_con_peak_infor,strcat(filename,'_con_MR'),'fig');
saveas(Cell_con_peak_infor,strcat(filename,'_con_MR'),'png');


