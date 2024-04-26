%%%%%%%% 有趋化作用的情况下，在趋化系数很小，该扩张模式和fisherwave类似时，计算侵染面积及各速度。
%%%%%%%% 首先计算纵向上Total_cell=0.9的位置FW_l为frontwave上特定的位置点；
%%%%%%%% 再在该位置处横向方向查找infected cell 的极值点为侵染区域的左右边界；
%%%%%%%% 拟合左右边界，从两边界的焦点开始至frontwave上FW_l范围内计算三角形面积，即为侵染面积；

clear all;
close all;


TimeT=220:250;
Diff_value =[2];%[60 90 120 150 180];

for i=1:length(Diff_value)
for m=1:length(TimeT)

    filename=strcat('initphage=1_GR=30_kap=2p5_PR=15_P0=0p05_Nk=5_Diff=',int2str(Diff_value(i)),'_',int2str(TimeT(m)));
    filename1=strcat('initphage=1_GR=30_kap=2.5_PR=15_P0=0.05_Nk=5_Diff=',int2str(Diff_value(i)),'_',int2str(TimeT(m)));

     load(strcat(filename1,'.mat'));
     %%  frontwave上infected cell半波长的计算
     Total_cell=Cell_den_S+Cell_den_I+Cell_den_R;
     infected_cell=Cell_den_I+Cell_den_R;
     
     ref_value=20;

  
     
     [~,FW_loca]=min(abs(Total_cell(ref_value:end,50)-1));  %计算fisherwave上total density为0.9的位置
      FW_l=FW_loca(end)+ref_value-1;
          
     %计算固定位置处侵染面积的左半波和右半波
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
%          
         
         
         
         
     
     
     [FWpeak_v,FWpeak_l]=findpeaks(infected_cell(FW_l,:));
     
     
      if ~isnan(length(FWpeak_l))
          
         [~,FW_Cell_Infect_half_loc_B]=min(abs(infected_cell(FW_l,1:FWpeak_l(1))-infect_ref_value));
         [~,FW_Cell_Infect_half_loc_E]=min(abs(infected_cell(FW_l,FWpeak_l(end):end)-infect_ref_value));
         
         FW_location_x(m)=FW_l;
         FW_location_y_left(m)=FW_Cell_Infect_half_loc_B;
         FW_location_y_right(m)=FW_Cell_Infect_half_loc_E+FWpeak_l(end)-1;

         FW_Infection_range(m)=(FW_location_y_right(m)-FW_location_y_left(m))*dx/10000; 
      
      end

     %%  Exp_MR 的计算  peak information of infected Cell
         fixed_exp_V=round(LY/dy)/2-50;
     
         [~,FW_exp_loca]=min(abs(Total_cell(ref_value:end,fixed_exp_V)-0.5));  %计算fisherwave上total density为0.9的位置
         FW_exp_l=FW_exp_loca(end)+ref_value-1;

         cell_infected_peak_location(m)=FW_exp_l*dx/10000; %cm
         cell_infected_peak_value(m)=Total_cell(FW_exp_l,fixed_exp_V);
     
     
   %% Con_MR的计算  peak information of  Mix Cell
     FW_con_l=FW_l;
     cell_con_peak_location(m)=FW_con_l*dx/10000; %cm
     cell_con_peak_value(m)=Total_cell(FW_con_l,100);
      
end
end

  %% migration velocity Cell 实验组
  real_time_h=TimeT;
cell_infected_peak_time=real_time_h;
cell_infected_peak_time(isnan(cell_infected_peak_location))=[];
cell_infected_peak_location(isnan(cell_infected_peak_location))=[];
MR_cell_infected_para=polyfit(cell_infected_peak_time,cell_infected_peak_location,1);
MR_cell_infected=MR_cell_infected_para(1)*cell_infected_peak_time+MR_cell_infected_para(2);
cell_infected_MR_V=MR_cell_infected_para(1);
% R2_Cell_infected_line(i,:)=norm(MR_cell_infected-mean(cell_infected_peak_location))^2/norm(cell_infected_peak_location-mean(cell_infected_peak_location))^2;
%% migration velocity Cell 对照组
cell_con_peak_time=real_time_h;
cell_con_peak_time(isnan(cell_con_peak_location))=[];
cell_con_peak_location(isnan(cell_con_peak_location))=[];
MR_cell_con_para=polyfit(cell_con_peak_time,cell_con_peak_location,1);
MR_cell_con=MR_cell_con_para(1)*cell_con_peak_time+MR_cell_con_para(2);
cell_con_MR_V=MR_cell_con_para(1);
% R2_Cell_con_line(i,:)=norm(MR_cell_con-mean(cell_con_peak_location))^2/norm(cell_con_peak_location-mean(cell_con_peak_location))^2;

%% infection speed  V//的计算
cell_infection_time=real_time_h;
Infection_speed_para=polyfit(cell_infection_time,FW_Infection_range,1);
Infection_speed_line=Infection_speed_para(1)*cell_infection_time+Infection_speed_para(2);
Infection_speed_V=Infection_speed_para(1);

%% 侵染边界拟合 左边界 横坐标是FW_location_y_left，纵坐标是FW_location_x（infected_area_left_line）
infected_area_left_para=polyfit(FW_location_y_left,FW_location_x,1);
infected_area_left_line=infected_area_left_para(1)*FW_location_y_left+infected_area_left_para(2);
%% 侵染边界拟合 右边界 横坐标是FW_location_y_right，纵坐标是FW_location_x（infected_area_right_line）
infected_area_right_para=polyfit(FW_location_y_right,FW_location_x,1);
infected_area_right_line=infected_area_right_para(1)*FW_location_y_right+infected_area_right_para(2);
%% 侵染面积的计算  已知纵坐标infection_X_range  求解横坐标left_location， 拟合曲线应反过来；横坐标的差值是每一行V-shape的面积差值
% Location_x=(infected_area_right_para(2)-infected_area_left_para(2))/(infected_area_left_para(1)-infected_area_right_para(1));
%  focal_point_temp=infected_area_left_para(1)*Location_x+infected_area_left_para(2);
%  focal_point=round(focal_point_temp);

infection_X_range=1:FW_l;

% infection_X_range=round(LXg0/dx):FW_l(length(FW_l));
for k=1:length(infection_X_range)
    left_location(k)=(infection_X_range(k)-infected_area_left_para(2))/infected_area_left_para(1);
    right_location(k)=(infection_X_range(k)-infected_area_right_para(2))/infected_area_right_para(1);
    infected_distance_X(k)=right_location(k)-left_location(k);
end  

[line_focus_v, line_focus_loc]=min(abs(round(left_location)-round(right_location)));

 if  ~isnan(line_focus_loc)
     infected_area=sum(infected_distance_X(line_focus_loc:end))*dx*dy;
     area_R=(length(infection_X_range)-line_focus_loc)*dx; 
     infected_ratio=2*infected_area/area_R/area_R;
 else
 infected_area=sum(infected_distance_X)*dx*dy;
 area_R=length(infection_X_range)*dx; 
 infected_ratio=infected_area/area_R/area_R;
 end 
 
 %% V-shape 角度的计算
 L_vect=[left_location(end)-left_location(1),infection_X_range(end)-infection_X_range(1)];
 R_vect=[right_location(end)-right_location(1),infection_X_range(end)-infection_X_range(1)];
 V_vect_angle=acos(dot(L_vect,R_vect)/(norm(L_vect)*norm( R_vect)));
 %% 图片 V-shape 边界及其拟合线
Cell_line=figure('position',[100 100 500 500]);
imagesc(infected_cell);hold on;
plot(FW_location_y_left,FW_location_x,'linewidth',4);hold on;
plot(FW_location_y_right,FW_location_x,'linewidth',4);
saveas(Cell_line,strcat(filename,'_Cell_line'),'fig');
saveas(Cell_line,strcat(filename,'_Cell_line'),'png');

 %% 图片 V-shape 边界及其拟合线
 Cell_simu_fit_line=figure('position',[100 100 500 500]);
imagesc(Total_cell);hold on;
 plot(FW_location_y_left,FW_location_x,'linewidth',4);hold on;
 plot(left_location,infection_X_range,'linewidth',4); hold on;
 plot(FW_location_y_right,FW_location_x,'linewidth',4);hold on;
 plot(right_location,infection_X_range,'linewidth',4); hold on;
 text(round(LY/dy)/2,round(LX/dx)/2,sprintf('infection-area-ratio=%.3f',infected_ratio),'FontSize',30);
 text(round(LY/dy)/3,round(LX/dx)/3,sprintf('IS-to-MR=%.3f',Infection_speed_V/cell_con_MR_V/2),'FontSize',30);
 text(round(LY/dy)/5,round(LX/dx)/5,sprintf('V-vect-angle=%.3f',V_vect_angle),'FontSize',30);
 saveas(Cell_simu_fit_line,strcat(filename,'_Cell_simu_fit_line'),'fig');
saveas(Cell_simu_fit_line,strcat(filename,'_Cell_simu_fit_line'),'png');

%% infection range 及其拟合线
FW_IR=figure('position',[100 100 500 500]);
plot(real_time_h,FW_Infection_range,'linewidth',4); hold on;
plot(cell_infection_time,Infection_speed_line,'r');hold on;
text(cell_infection_time(end-10),FW_Infection_range(end-10),sprintf('IS=%.3f',Infection_speed_V),'FontSize',30);

xlabel('time');
ylabel('infection range');
% legend('chi=200');
title('front wave infection range');
set(gca,'linewidth',3,'FontSize',30,'LineWidth',3);
saveas(FW_IR,strcat(filename,'_Infection_speed'),'fig');
saveas(FW_IR,strcat(filename,'_Infection_speed'),'png');
%%  实验组 速度及其拟合
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

%% 对照组速度及其拟合
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

% %%
% C_IR=figure('position',[100 100 1500 900]);
% plot(TimeT,C_Infection_range,'linewidth',4); 
% xlabel('time');
% ylabel('infection range');
% title('Cell_den(400,:) infection range');
% legend('chi=60','chi=90','chi=120','chi=150','chi=180')
% set(gca,'linewidth',3,'FontSize',30,'LineWidth',3);
% saveas(C_IR,strcat(filename,'_C_IR'),'fig');
% saveas(C_IR,strcat(filename,'_C_IR'),'png');

%%
% FW_Cell_1_IR=figure('position',[100 100 1500 900]);
% 
% % h1=plot(Cell_den(FW_l(end),:,1)+Cell_den(FW_l(end),:,2),'linewidth',4,'Color','black');hold on;
% % h2=plot(Cell_den(FW_l(end),:,1),'linewidth',3,'Color',[0 0.5 0]);hold on;
% h3=plot(Total_cell(FW_l(end),:),'linewidth',3,'Color','red');hold on;
% plot(FW_Cell_1_half_loc_B,Total_cell(FW_l(end),FW_Cell_1_half_loc_B),'*','Color','black');hold on;
%  plot(FW_Cell_1_half_loc_E,Total_cell(FW_l(end),FW_Cell_1_half_loc_E),'*','Color','black');hold on;
% % xlim([0 320]);
% % ylim([0 1.3]);
% xlabel('Distance');
% % ylabel('OD_6_0_0');
% % set(gca,'xtick',[0 100 200 300]);
% % set(gca,'ytick',[0 0.4 0.8 1.2]);
% set(gca,'linewidth',3,'FontSize',30,'LineWidth',3);
% % l1=legend([h1 h2 h3],{'total cell','uninfected cell','infected cell'},'location','Northwest');
% % set(l1,'Fontname', 'Times New Roman','FontWeight','bold','FontSize',30,'box','off');
% % text(250,max(Nutr(400,300))+0.04,'rel.nutrient','FontSize',30);
% title('Front wave infection range')
% saveas(FW_Cell_1_IR,strcat(filename,'_cell'),'fig');
% saveas(FW_Cell_1_IR,strcat(filename,'_cell'),'png');

%  %%
% h=figure('position',[100 100 1500 900]);
% c=area(Nutr(400,:),'LineStyle',':');c.FaceColor = [1 0.7 0.2];hold on;
% h1=plot(Cell_den(400,:,1)+Cell_den(400,:,2),'linewidth',4,'Color','black');hold on;
% h2=plot(Cell_den(400,:,1),'linewidth',3,'Color',[0 0.5 0]);hold on;
% h3=plot(Cell_den(400,:,2),'linewidth',3,'Color','red');hold on;
% plot(Cell_2_half_loc_B,Cell_den(400,Cell_2_half_loc_B,2),'*','Color','black');hold on;
% plot(Cell_2_half_loc_E,Cell_den(400,Cell_2_half_loc_E,2),'*','Color','black');hold on;
% xlim([0 320]);
% ylim([0 1.3]);
% xlabel('Distance');
% % ylabel('OD_6_0_0');
% set(gca,'xtick',[0 100 200 300]);
% set(gca,'ytick',[0 0.4 0.8 1.2]);
% set(gca,'linewidth',3,'FontSize',30,'LineWidth',3);
% l1=legend([h1 h2 h3],{'total cell','uninfected cell','infected cell'},'location','Northwest');
% set(l1,'Fontname', 'Times New Roman','FontWeight','bold','FontSize',30,'box','off');
% text(250,max(Nutr(400,300))+0.04,'rel.nutrient','FontSize',30);
% title('Cell_den(400,:) infection range')
% saveas(h,strcat(filename,'_cell'),'fig');
% saveas(h,strcat(filename,'_cell'),'png');
