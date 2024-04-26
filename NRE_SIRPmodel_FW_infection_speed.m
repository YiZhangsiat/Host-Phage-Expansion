%%%%%%%% 有趋化作用的情况下，在趋化系数很小，该扩张模式和fisherwave类似时，计算侵染面积及各速度。
%%%%%%%% 首先计算纵向上Total_cell=0.9的位置FW_l为frontwave上特定的位置点；
%%%%%%%% 再在该位置处横向方向查找infected cell 的极值点为侵染区域的左右边界；
%%%%%%%% 拟合左右边界，从两边界的焦点开始至frontwave上FW_l范围内计算三角形面积，即为侵染面积；

clear all;
close all;


TimeT=20:40;
eta_value =[8];%[60 90 120 150 180];

for i=1:length(eta_value)
for m=1:length(TimeT)
 
    filename=strcat('C_tar35Simu_initP=1E-3_diff=30_K=3p5_PR=10_chi=240_eta=',int2str(eta_value(i)),'_',int2str(TimeT(m)));
    filename1=strcat('C_tar35Simu_initP=1E-3_diff=30_K=3.5_PR=10_chi=240_eta=',int2str(eta_value(i)),'_',int2str(TimeT(m)));

     load(strcat(filename1,'.mat'));
     %%  frontwave上infected cell半波长的计算
     Total_cell=Cell_den_S+Cell_den_I+Cell_den_R;
     infected_cell=Cell_den_I+Cell_den_R;
     

     [Nx,Ny]=size(Total_cell);

     Nx0=round(Nx/2); 
     Ny0=round(Ny/2);
     

     %%  Exp_MR 的计算  peak information of infected Cell
       [cell_infected_maxv,cell_infected_maxl]=findpeaks(Total_cell(Nx0:end,Ny0));  % 计算实验组front wave所处的位置 
  if  ~isnan(cell_infected_maxv)
  cell_infected_peak_value(m)=cell_infected_maxv(end);
  FW_exp_l=cell_infected_maxl(end)+Nx0-1;
  cell_infected_peak_location(m)=FW_exp_l*dx/10000; %cm
  else
  cell_infected_peak_value(m)=nan;
  FW_exp_l=nan; %cm
  cell_infected_peak_location(m)=nan; %cm
  end
     
     
 
     
     
   %% Con_MR的计算  peak information of  Mix Cell
   
       [cell_con_maxv,cell_con_maxl]=findpeaks(Total_cell(Nx0,Ny0:end));  
   
      FW_con_l=cell_con_maxl(end)+Ny0-1;

     cell_con_peak_location(m)=FW_con_l*dx/10000; %cm
     cell_con_peak_value(m)=Total_cell(Nx0,FW_con_l);
      
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


%%  实验组 速度及其拟合
Cell_Exp_peak_infor=figure('position',[100 100 500 500]);
subplot(2,1,1);plot(real_time_h,cell_infected_peak_value); 
text(real_time_h(5),cell_infected_peak_value(5),sprintf('peak value=%.3f',max(cell_infected_peak_value)),'FontSize',30);
% legend('Cell 1 peak value');
xlabel('time');
ylabel('Cell infected peak value');
title('cell infected front peak')
subplot(2,1,2);plot(cell_infected_peak_time,cell_infected_peak_location); 
hold on;plot(cell_infected_peak_time,MR_cell_infected,'r');
text(real_time_h(5),cell_infected_peak_location(5),sprintf('MR=%.3f',cell_infected_MR_V),'FontSize',30);
% legend('Cell 1 major peak location','Cell 1 line');
xlabel('time');
ylabel('Cell infected peak location');

saveas(Cell_Exp_peak_infor,strcat(filename,'_Exp_MR'),'fig');
saveas(Cell_Exp_peak_infor,strcat(filename,'_Exp_MR'),'png');

%% 对照组速度及其拟合
Cell_con_peak_infor=figure('position',[100 100 500 500]);
subplot(2,1,1);plot(real_time_h,cell_con_peak_value); 
text(real_time_h(5),cell_con_peak_value(5),sprintf('peak value=%.3f',max(cell_con_peak_value)),'FontSize',30);
% legend('Cell 1 peak value');
xlabel('time');
ylabel('cell con peak_value');
title('con cell front peak')
subplot(2,1,2);plot(cell_con_peak_time,cell_con_peak_location); 
hold on;plot(cell_con_peak_time,MR_cell_con,'r');
text(real_time_h(5),cell_con_peak_location(end),sprintf('MR=%.3f',cell_con_MR_V),'FontSize',30);
% legend('Cell 1 major peak location','Cell 1 line');
xlabel('time');
ylabel('Cell con peak location');
saveas(Cell_con_peak_infor,strcat(filename,'_con_MR'),'fig');
saveas(Cell_con_peak_infor,strcat(filename,'_con_MR'),'png');



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


