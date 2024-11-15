 clear all;
close all;

TimeT=[21];

PR=[400];
 filename=strcat('newestmodel_constInfection_initphage=1E-1_Dp=0_attr=100_chi=300_phagePR=',int2str(PR),'_',int2str(TimeT));
 filename1=strcat('newestmodel_constInfection_initphage=1E-1_Dp=0_attr=100_chi=300_phagePR=',int2str(PR),'_',int2str(TimeT));
load(strcat(filename1,'.mat'));

     %%  
  cal_region=round(80000/250)/2;
 Total_cell1=Cell_den_S+Cell_den_I+Cell_den_R;
 Total_cell=Cell_den_S+Cell_den_I+Cell_den_R;
 Infected_Cell=Cell_den_I+Cell_den_R;
 
 
 [Nx,Ny]=size(Total_cell);

Nx0=round(Nx/2); 
Ny0=round(Ny/2);

for i=1:Nx
    for j=1:Ny
        if (i-Nx0)^2+(j-Ny0)^2>cal_region^2
      Total_cell(i,j)=0;
      Infected_Cell(i,j)=0;
      Phag(i,j)=0;
        end
    end
end

Infected_Cell(Infected_Cell<0.1)=0;
Cell_Area_logical_temp=(Infected_Cell~=0);
Cell_Area_total_pixel_1=sum(Cell_Area_logical_temp'); 
line_cell_area=sum(Cell_Area_total_pixel_1~=0);    
Cell_Area_total_pixel=sum(Cell_Area_total_pixel_1)-line_cell_area;  

Phag(Phag<0.05)=0;
Phage_Area_logical_temp=(Phag~=0);
Phage_Area_total_pixel_1=sum(Phage_Area_logical_temp'); 
line_Phage_area=sum(Phage_Area_total_pixel_1~=0);    
Phage_Area_total_pixel=sum(Phage_Area_total_pixel_1)-line_Phage_area;  
unit_pixel_area=dx*dy/10000/10000;      
Cell_Area=Cell_Area_total_pixel*unit_pixel_area;  
Phage_Area=Phage_Area_total_pixel*unit_pixel_area;  

angle_phage=2*Phage_Area/3/3;
angle_cell=2*Cell_Area/3/3;

%%
Fig1=figure('position',[100 100 550 550]);
imagesc(Phag);axis square;
text(Nx0,Ny0,sprintf('infected phage area=%.3f',max(Phage_Area)),'Color','red','FontSize',18);
text(1,400,sprintf('angle=%.3f',max(angle_phage)),'Color','red','FontSize',18);
title('phage');
saveas(Fig1,strcat(filename,'_infected_Phage_Area'),'png');
