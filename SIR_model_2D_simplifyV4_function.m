% clear all;
% close all;
function []=SIR_model_2D_simplifyV4_function(name,eta)
%% define parameters
dt=1;                  
dx=100;                 
dy=100;
Time=60*60*40;        
T=round(Time/dt);       
LX=120000;            
LY=120000;
LX0=LX/2; 
LY0=LY/2;
Rb1=2000;
Rb2=1000;

LXg0=LX0+10000; %50000;
LYg0=LY/2;
Rg1=2000; 
Rg2=1000;


K1=3.5;               
K2=1000;              
Sk=1;                
Nk=5;                
Nk_2=Nk*Nk;
Yn=0.06;


Dn=800;               
Da=800;              
G0=4/60;        
Nutr_init=30;            
Attr_init=100;       
bacter_init=0.2;
phage_init=0.001;

Diff=30;
chi=180;    
lamnp0=log(2)/30/60;
kap0=2.5;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
P0=0.05;

% eta=6;
PT0=P0*P0;
beta1=0;
beta2=0.8;
theta=5E-4;

%% initial conditions




Cell_den_S=zeros(round(LX/dx),round(LY/dy));
Cell_den_I=zeros(round(LX/dx),round(LY/dy));
Cell_den_R=zeros(round(LX/dx),round(LY/dy));

Nutr=ones(round(LX/dx),round(LY/dy))*Nutr_init;
Attr=ones(round(LX/dx),round(LY/dy))*Attr_init;
Phag=zeros(round(LX/dx),round(LY/dy));

tic

% Cell_den_S(round(LX0/dx)-2:round(LX0/dx)+2,:)=bacter_init;
for x=1:round(LX/dx)
    for y=1:round(LY/dy)
        RR1=Rb1^2;
        RR2=Rb2^2;
        r2=((x-1)*dx-LX0)^2+((y-1)*dy-LY0)^2;
        if r2<=RR1
            bb(x,y)=bacter_init*exp(-r2/RR2);
        else
            bb(x,y)=0;
        end
    end
end  
Cell_den_S=bb;




for x=1:round(LX/dx)
    for y=1:round(LY/dy)
        RRg1=Rg1^2;
        RRg2=Rg2^2;
        r2=((x-1)*dx-LXg0)^2+((y-1)*dy-LYg0)^2;
        if r2<=RRg1
            pp(x,y)=phage_init*exp(-r2/RRg2);
        else
            pp(x,y)=0;
        end
    end
end  
Phag=pp;

%%
Cell_den_S_temp_a=Cell_den_S;Cell_den_S_temp_b=Cell_den_S;Cell_den_S_temp_c=Cell_den_S;Cell_den_S_temp_d=Cell_den_S;
Cell_den_I_temp_a=Cell_den_I;Cell_den_I_temp_b=Cell_den_I;Cell_den_I_temp_c=Cell_den_I;Cell_den_I_temp_d=Cell_den_I;
Cell_den_R_temp_a=Cell_den_R;Cell_den_R_temp_b=Cell_den_R;Cell_den_R_temp_c=Cell_den_R;Cell_den_R_temp_d=Cell_den_R;

Attr_temp_a=Attr;Attr_temp_b=Attr;Attr_temp_c=Attr;Attr_temp_d=Attr;
Nutr_temp_a=Nutr;Nutr_temp_b=Nutr;Nutr_temp_c=Nutr;Nutr_temp_d=Nutr;
Phag_temp_a=Phag;Phag_temp_b=Phag;Phag_temp_c=Phag;Phag_temp_d=Phag;
g1_S=Cell_den_S;g2_S=Cell_den_S;
g1_I=Cell_den_I;g2_I=Cell_den_I;
g1_R=Cell_den_R;g2_R=Cell_den_R;
%%
for t=1:T

    
    Cell_den_S_temp=Cell_den_S;
    Cell_den_S_temp(Cell_den_S_temp<=1E-10)=0;
    Cell_den_I_temp=Cell_den_I;
    Cell_den_I_temp(Cell_den_I_temp<=1E-10)=0;
    Cell_den_R_temp=Cell_den_R;
    Cell_den_R_temp(Cell_den_R_temp<=1E-10)=0; 
    Attr_temp=Attr;
    Attr_temp(Attr_temp<=0)=0;
    Nutr_temp=Nutr; 
    Nutr_temp(Nutr_temp<=0)=0;
    Phag_temp=Phag; 
    Phag_temp(Phag_temp<=1E-10)=0;
    
    
    
    lambda_np_temp=lamnp0./(1+Nk_2./Nutr_temp./Nutr_temp);
    kappa_temp=kap0./(1+PT0./Phag_temp./Phag_temp);
    Gamma_temp=G0./(1+Sk./Attr_temp);
    
    S_growth_temp=lambda_np_temp.*Cell_den_S_temp-kappa_temp.*lambda_np_temp.*Cell_den_S_temp;
    I_growth_temp=beta1*lambda_np_temp.*Cell_den_I_temp+kappa_temp.*lambda_np_temp.*Cell_den_S_temp-theta*Cell_den_I_temp;
    R_growth_temp=beta2*lambda_np_temp.*Cell_den_R_temp+theta*Cell_den_I_temp;
    
    nutr_consump_temp=lambda_np_temp.*(Cell_den_S_temp+Cell_den_I_temp+Cell_den_R_temp)./Yn;
    attr_consump_temp=Gamma_temp.*(Cell_den_S_temp+Cell_den_I_temp+Cell_den_R_temp);
    phage_product_temp=eta*((1-beta1)*lambda_np_temp.*Cell_den_I_temp+(1-beta2)*lambda_np_temp.*Cell_den_R_temp);
    
    
    f=chi.*log((1+ Attr_temp/K1)./(1+ Attr_temp/K2));% free energy
    g1_S(2:end-1,:)=Cell_den_S_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_S(:,2:end-1)=Cell_den_S_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_S(1,:)=-g1_S(2,:);g1_S(end,:)=-g1_S(end-1,:);
    g2_S(:,1)=-g2_S(:,2);g2_S(:,end)=-g2_S(:,end-1); 
    
    g1_I(2:end-1,:)=Cell_den_I_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_I(:,2:end-1)=Cell_den_I_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_I(1,:)=-g1_I(2,:);g1_I(end,:)=-g1_I(end-1,:);
    g2_I(:,1)=-g2_I(:,2);g2_I(:,end)=-g2_I(:,end-1); 
    
    g1_R(2:end-1,:)=Cell_den_R_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_R(:,2:end-1)=Cell_den_R_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_R(1,:)=-g1_R(2,:);g1_R(end,:)=-g1_R(end-1,:);
    g2_R(:,1)=-g2_R(:,2);g2_R(:,end)=-g2_R(:,end-1); 


    
    Cell_den_S_temp_a(2:end-1,2:end-1)=Diff.*((Cell_den_S_temp(1:end-2,2:end-1)+Cell_den_S_temp(3:end,2:end-1)-2*Cell_den_S_temp(2:end-1,2:end-1))/dx^2+(Cell_den_S_temp(2:end-1,1:end-2)+Cell_den_S_temp(2:end-1,3:end)-2*Cell_den_S_temp(2:end-1,2:end-1))/dy^2)-(g1_S(3:end,2:end-1)-g1_S(1:end-2,2:end-1))/2/dx-(g2_S(2:end-1,3:end)-g2_S(2:end-1,1:end-2))/2/dy+S_growth_temp(2:end-1,2:end-1);
    Cell_den_I_temp_a(2:end-1,2:end-1)=Diff.*((Cell_den_I_temp(1:end-2,2:end-1)+Cell_den_I_temp(3:end,2:end-1)-2*Cell_den_I_temp(2:end-1,2:end-1))/dx^2+(Cell_den_I_temp(2:end-1,1:end-2)+Cell_den_I_temp(2:end-1,3:end)-2*Cell_den_I_temp(2:end-1,2:end-1))/dy^2)-(g1_I(3:end,2:end-1)-g1_I(1:end-2,2:end-1))/2/dx-(g2_I(2:end-1,3:end)-g2_I(2:end-1,1:end-2))/2/dy+I_growth_temp(2:end-1,2:end-1);
    Cell_den_R_temp_a(2:end-1,2:end-1)=Diff.*((Cell_den_R_temp(1:end-2,2:end-1)+Cell_den_R_temp(3:end,2:end-1)-2*Cell_den_R_temp(2:end-1,2:end-1))/dx^2+(Cell_den_R_temp(2:end-1,1:end-2)+Cell_den_R_temp(2:end-1,3:end)-2*Cell_den_R_temp(2:end-1,2:end-1))/dy^2)-(g1_R(3:end,2:end-1)-g1_R(1:end-2,2:end-1))/2/dx-(g2_R(2:end-1,3:end)-g2_R(2:end-1,1:end-2))/2/dy+R_growth_temp(2:end-1,2:end-1);
    
    Attr_temp_a(2:end-1,2:end-1)=Da*((Attr_temp(1:end-2,2:end-1)+Attr_temp(3:end,2:end-1)-2*Attr_temp(2:end-1,2:end-1))/dx^2+(Attr_temp(2:end-1,1:end-2)+Attr_temp(2:end-1,3:end)-2*Attr_temp(2:end-1,2:end-1))/dy^2)-attr_consump_temp(2:end-1,2:end-1);
    Nutr_temp_a(2:end-1,2:end-1)=Dn*((Nutr_temp(1:end-2,2:end-1)+Nutr_temp(3:end,2:end-1)-2*Nutr_temp(2:end-1,2:end-1))/dx^2+(Nutr_temp(2:end-1,1:end-2)+Nutr_temp(2:end-1,3:end)-2*Nutr_temp(2:end-1,2:end-1))/dy^2)-nutr_consump_temp(2:end-1,2:end-1);
    Phag_temp_a(2:end-1,2:end-1)=phage_product_temp(2:end-1,2:end-1);
    
    
 
    Cell_den_S_temp_a(1,:)= Cell_den_S_temp_a(2,:);  Cell_den_S_temp_a(:,1)= Cell_den_S_temp_a(:,2); 
    Cell_den_S_temp_a(end,:)= Cell_den_S_temp_a(end-1,:);   Cell_den_S_temp_a(:,end)= Cell_den_S_temp_a(:,end-1); 
    Cell_den_I_temp_a(1,:)= Cell_den_I_temp_a(2,:);  Cell_den_I_temp_a(:,1)= Cell_den_I_temp_a(:,2); 
    Cell_den_I_temp_a(end,:)= Cell_den_I_temp_a(end-1,:);   Cell_den_I_temp_a(:,end)= Cell_den_I_temp_a(:,end-1); 
    Cell_den_R_temp_a(1,:)= Cell_den_R_temp_a(2,:);  Cell_den_R_temp_a(:,1)= Cell_den_R_temp_a(:,2); 
    Cell_den_R_temp_a(end,:)= Cell_den_R_temp_a(end-1,:);   Cell_den_R_temp_a(:,end)= Cell_den_R_temp_a(:,end-1); 
    
    Attr_temp_a(1,:)=Attr_temp_a(2,:); Attr_temp_a(:,1)=Attr_temp_a(:,2); 
    Attr_temp_a(end,:)=Attr_temp_a(end-1,:);  Attr_temp_a(:,end)=Attr_temp_a(:,end-1); 
    Nutr_temp_a(1,:)=Nutr_temp_a(2,:); Nutr_temp_a(:,1)=Nutr_temp_a(:,2); 
    Nutr_temp_a(end,:)=Nutr_temp_a(end-1,:);  Nutr_temp_a(:,end)=Nutr_temp_a(:,end-1); 
    Phag_temp_a(1,:)=Phag_temp_a(2,:); Phag_temp_a(:,1)=Phag_temp_a(:,2); 
    Phag_temp_a(end,:)=Phag_temp_a(end-1,:);  Phag_temp_a(:,end)=Phag_temp_a(:,end-1);
    
  
 
    %%  
    Cell_den_S_temp=Cell_den_S+Cell_den_S_temp_a*dt/2;
    Cell_den_S_temp(Cell_den_S_temp<=1E-10)=0;
    Cell_den_I_temp=Cell_den_I+Cell_den_I_temp_a*dt/2;
    Cell_den_I_temp(Cell_den_I_temp<=1E-10)=0;
    Cell_den_R_temp=Cell_den_R+Cell_den_R_temp_a*dt/2;
    Cell_den_R_temp(Cell_den_R_temp<=1E-10)=0;  
    Attr_temp=Attr+Attr_temp_a*dt/2;
    Attr_temp(Attr_temp<=0)=0;
    Nutr_temp=Nutr+Nutr_temp_a*dt/2; 
    Nutr_temp(Nutr_temp<=0)=0;
    Phag_temp=Phag+Phag_temp_a*dt/2; 
    Phag_temp(Phag_temp<=1E-10)=0;
    
  
    lambda_np_temp=lamnp0./(1+Nk_2./Nutr_temp./Nutr_temp);
    kappa_temp=kap0./(1+PT0./Phag_temp./Phag_temp);
    Gamma_temp=G0./(1+Sk./Attr_temp);
    
    S_growth_temp=lambda_np_temp.*Cell_den_S_temp-kappa_temp.*lambda_np_temp.*Cell_den_S_temp;
    I_growth_temp=beta1*lambda_np_temp.*Cell_den_I_temp+kappa_temp.*lambda_np_temp.*Cell_den_S_temp-theta*Cell_den_I_temp;
    R_growth_temp=beta2*lambda_np_temp.*Cell_den_R_temp+theta*Cell_den_I_temp;
    
    nutr_consump_temp=lambda_np_temp.*(Cell_den_S_temp+Cell_den_I_temp+Cell_den_R_temp)./Yn;
    attr_consump_temp=Gamma_temp.*(Cell_den_S_temp+Cell_den_I_temp+Cell_den_R_temp);
    phage_product_temp=eta*((1-beta1)*lambda_np_temp.*Cell_den_I_temp+(1-beta2)*lambda_np_temp.*Cell_den_R_temp);
    
    
    f=chi.*log((1+ Attr_temp/K1)./(1+ Attr_temp/K2));% free energy
    g1_S(2:end-1,:)=Cell_den_S_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_S(:,2:end-1)=Cell_den_S_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_S(1,:)=-g1_S(2,:);g1_S(end,:)=-g1_S(end-1,:);
    g2_S(:,1)=-g2_S(:,2);g2_S(:,end)=-g2_S(:,end-1); 
    
    g1_I(2:end-1,:)=Cell_den_I_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_I(:,2:end-1)=Cell_den_I_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_I(1,:)=-g1_I(2,:);g1_I(end,:)=-g1_I(end-1,:);
    g2_I(:,1)=-g2_I(:,2);g2_I(:,end)=-g2_I(:,end-1); 
    
    g1_R(2:end-1,:)=Cell_den_R_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_R(:,2:end-1)=Cell_den_R_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_R(1,:)=-g1_R(2,:);g1_R(end,:)=-g1_R(end-1,:);
    g2_R(:,1)=-g2_R(:,2);g2_R(:,end)=-g2_R(:,end-1); 


    
    Cell_den_S_temp_b(2:end-1,2:end-1)=Diff.*((Cell_den_S_temp(1:end-2,2:end-1)+Cell_den_S_temp(3:end,2:end-1)-2*Cell_den_S_temp(2:end-1,2:end-1))/dx^2+(Cell_den_S_temp(2:end-1,1:end-2)+Cell_den_S_temp(2:end-1,3:end)-2*Cell_den_S_temp(2:end-1,2:end-1))/dy^2)-(g1_S(3:end,2:end-1)-g1_S(1:end-2,2:end-1))/2/dx-(g2_S(2:end-1,3:end)-g2_S(2:end-1,1:end-2))/2/dy+S_growth_temp(2:end-1,2:end-1);
    Cell_den_I_temp_b(2:end-1,2:end-1)=Diff.*((Cell_den_I_temp(1:end-2,2:end-1)+Cell_den_I_temp(3:end,2:end-1)-2*Cell_den_I_temp(2:end-1,2:end-1))/dx^2+(Cell_den_I_temp(2:end-1,1:end-2)+Cell_den_I_temp(2:end-1,3:end)-2*Cell_den_I_temp(2:end-1,2:end-1))/dy^2)-(g1_I(3:end,2:end-1)-g1_I(1:end-2,2:end-1))/2/dx-(g2_I(2:end-1,3:end)-g2_I(2:end-1,1:end-2))/2/dy+I_growth_temp(2:end-1,2:end-1);
    Cell_den_R_temp_b(2:end-1,2:end-1)=Diff.*((Cell_den_R_temp(1:end-2,2:end-1)+Cell_den_R_temp(3:end,2:end-1)-2*Cell_den_R_temp(2:end-1,2:end-1))/dx^2+(Cell_den_R_temp(2:end-1,1:end-2)+Cell_den_R_temp(2:end-1,3:end)-2*Cell_den_R_temp(2:end-1,2:end-1))/dy^2)-(g1_R(3:end,2:end-1)-g1_R(1:end-2,2:end-1))/2/dx-(g2_R(2:end-1,3:end)-g2_R(2:end-1,1:end-2))/2/dy+R_growth_temp(2:end-1,2:end-1);
    
    Attr_temp_b(2:end-1,2:end-1)=Da*((Attr_temp(1:end-2,2:end-1)+Attr_temp(3:end,2:end-1)-2*Attr_temp(2:end-1,2:end-1))/dx^2+(Attr_temp(2:end-1,1:end-2)+Attr_temp(2:end-1,3:end)-2*Attr_temp(2:end-1,2:end-1))/dy^2)-attr_consump_temp(2:end-1,2:end-1);
    Nutr_temp_b(2:end-1,2:end-1)=Dn*((Nutr_temp(1:end-2,2:end-1)+Nutr_temp(3:end,2:end-1)-2*Nutr_temp(2:end-1,2:end-1))/dx^2+(Nutr_temp(2:end-1,1:end-2)+Nutr_temp(2:end-1,3:end)-2*Nutr_temp(2:end-1,2:end-1))/dy^2)-nutr_consump_temp(2:end-1,2:end-1);
    Phag_temp_b(2:end-1,2:end-1)=phage_product_temp(2:end-1,2:end-1);
    
    
    
 
    Cell_den_S_temp_b(1,:)= Cell_den_S_temp_b(2,:);  Cell_den_S_temp_b(:,1)= Cell_den_S_temp_b(:,2); 
    Cell_den_S_temp_b(end,:)= Cell_den_S_temp_b(end-1,:);   Cell_den_S_temp_b(:,end)= Cell_den_S_temp_b(:,end-1); 
    Cell_den_I_temp_b(1,:)= Cell_den_I_temp_b(2,:);  Cell_den_I_temp_b(:,1)= Cell_den_I_temp_b(:,2); 
    Cell_den_I_temp_b(end,:)= Cell_den_I_temp_b(end-1,:);   Cell_den_I_temp_b(:,end)= Cell_den_I_temp_b(:,end-1); 
    Cell_den_R_temp_b(1,:)= Cell_den_R_temp_b(2,:);  Cell_den_R_temp_b(:,1)= Cell_den_R_temp_b(:,2); 
    Cell_den_R_temp_b(end,:)= Cell_den_R_temp_b(end-1,:);   Cell_den_R_temp_b(:,end)= Cell_den_R_temp_b(:,end-1); 
        
    Attr_temp_b(1,:)=Attr_temp_b(2,:); Attr_temp_b(:,1)=Attr_temp_b(:,2); 
    Attr_temp_b(end,:)=Attr_temp_b(end-1,:);  Attr_temp_b(:,end)=Attr_temp_b(:,end-1); 
    Nutr_temp_b(1,:)=Nutr_temp_b(2,:); Nutr_temp_b(:,1)=Nutr_temp_b(:,2); 
    Nutr_temp_b(end,:)=Nutr_temp_b(end-1,:);  Nutr_temp_b(:,end)=Nutr_temp_b(:,end-1); 
    Phag_temp_b(1,:)=Phag_temp_b(2,:); Phag_temp_b(:,1)=Phag_temp_b(:,2); 
    Phag_temp_b(end,:)=Phag_temp_b(end-1,:);  Phag_temp_b(:,end)=Phag_temp_b(:,end-1);   
    
    
    %%
    Cell_den_S_temp=Cell_den_S+Cell_den_S_temp_b*dt/2;
    Cell_den_S_temp(Cell_den_S_temp<=1E-10)=0;
    Cell_den_I_temp=Cell_den_I+Cell_den_I_temp_b*dt/2;
    Cell_den_I_temp(Cell_den_I_temp<=1E-10)=0;
    Cell_den_R_temp=Cell_den_R+Cell_den_R_temp_b*dt/2;
    Cell_den_R_temp(Cell_den_R_temp<=1E-10)=0;   
    Attr_temp=Attr+Attr_temp_b*dt/2;
    Attr_temp(Attr_temp<=0)=0;
    Nutr_temp=Nutr+Nutr_temp_b*dt/2; 
    Nutr_temp(Nutr_temp<=0)=0;
    Phag_temp=Phag+Phag_temp_b*dt/2; 
    Phag_temp(Phag_temp<=1E-10)=0;
  
    
    lambda_np_temp=lamnp0./(1+Nk_2./Nutr_temp./Nutr_temp);
    kappa_temp=kap0./(1+PT0./Phag_temp./Phag_temp);
    Gamma_temp=G0./(1+Sk./Attr_temp);
    
    S_growth_temp=lambda_np_temp.*Cell_den_S_temp-kappa_temp.*lambda_np_temp.*Cell_den_S_temp;
    I_growth_temp=beta1*lambda_np_temp.*Cell_den_I_temp+kappa_temp.*lambda_np_temp.*Cell_den_S_temp-theta*Cell_den_I_temp;
    R_growth_temp=beta2*lambda_np_temp.*Cell_den_R_temp+theta*Cell_den_I_temp;
    
    nutr_consump_temp=lambda_np_temp.*(Cell_den_S_temp+Cell_den_I_temp+Cell_den_R_temp)./Yn;
    attr_consump_temp=Gamma_temp.*(Cell_den_S_temp+Cell_den_I_temp+Cell_den_R_temp);
    phage_product_temp=eta*((1-beta1)*lambda_np_temp.*Cell_den_I_temp+(1-beta2)*lambda_np_temp.*Cell_den_R_temp);
    
    
    f=chi.*log((1+ Attr_temp/K1)./(1+ Attr_temp/K2));% free energy
    g1_S(2:end-1,:)=Cell_den_S_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_S(:,2:end-1)=Cell_den_S_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_S(1,:)=-g1_S(2,:);g1_S(end,:)=-g1_S(end-1,:);
    g2_S(:,1)=-g2_S(:,2);g2_S(:,end)=-g2_S(:,end-1); 
    
    g1_I(2:end-1,:)=Cell_den_I_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_I(:,2:end-1)=Cell_den_I_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_I(1,:)=-g1_I(2,:);g1_I(end,:)=-g1_I(end-1,:);
    g2_I(:,1)=-g2_I(:,2);g2_I(:,end)=-g2_I(:,end-1); 
    
    g1_R(2:end-1,:)=Cell_den_R_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_R(:,2:end-1)=Cell_den_R_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_R(1,:)=-g1_R(2,:);g1_R(end,:)=-g1_R(end-1,:);
    g2_R(:,1)=-g2_R(:,2);g2_R(:,end)=-g2_R(:,end-1); 


    
    Cell_den_S_temp_c(2:end-1,2:end-1)=Diff.*((Cell_den_S_temp(1:end-2,2:end-1)+Cell_den_S_temp(3:end,2:end-1)-2*Cell_den_S_temp(2:end-1,2:end-1))/dx^2+(Cell_den_S_temp(2:end-1,1:end-2)+Cell_den_S_temp(2:end-1,3:end)-2*Cell_den_S_temp(2:end-1,2:end-1))/dy^2)-(g1_S(3:end,2:end-1)-g1_S(1:end-2,2:end-1))/2/dx-(g2_S(2:end-1,3:end)-g2_S(2:end-1,1:end-2))/2/dy+S_growth_temp(2:end-1,2:end-1);
    Cell_den_I_temp_c(2:end-1,2:end-1)=Diff.*((Cell_den_I_temp(1:end-2,2:end-1)+Cell_den_I_temp(3:end,2:end-1)-2*Cell_den_I_temp(2:end-1,2:end-1))/dx^2+(Cell_den_I_temp(2:end-1,1:end-2)+Cell_den_I_temp(2:end-1,3:end)-2*Cell_den_I_temp(2:end-1,2:end-1))/dy^2)-(g1_I(3:end,2:end-1)-g1_I(1:end-2,2:end-1))/2/dx-(g2_I(2:end-1,3:end)-g2_I(2:end-1,1:end-2))/2/dy+I_growth_temp(2:end-1,2:end-1);
    Cell_den_R_temp_c(2:end-1,2:end-1)=Diff.*((Cell_den_R_temp(1:end-2,2:end-1)+Cell_den_R_temp(3:end,2:end-1)-2*Cell_den_R_temp(2:end-1,2:end-1))/dx^2+(Cell_den_R_temp(2:end-1,1:end-2)+Cell_den_R_temp(2:end-1,3:end)-2*Cell_den_R_temp(2:end-1,2:end-1))/dy^2)-(g1_R(3:end,2:end-1)-g1_R(1:end-2,2:end-1))/2/dx-(g2_R(2:end-1,3:end)-g2_R(2:end-1,1:end-2))/2/dy+R_growth_temp(2:end-1,2:end-1);
    
    Attr_temp_c(2:end-1,2:end-1)=Da*((Attr_temp(1:end-2,2:end-1)+Attr_temp(3:end,2:end-1)-2*Attr_temp(2:end-1,2:end-1))/dx^2+(Attr_temp(2:end-1,1:end-2)+Attr_temp(2:end-1,3:end)-2*Attr_temp(2:end-1,2:end-1))/dy^2)-attr_consump_temp(2:end-1,2:end-1);
    Nutr_temp_c(2:end-1,2:end-1)=Dn*((Nutr_temp(1:end-2,2:end-1)+Nutr_temp(3:end,2:end-1)-2*Nutr_temp(2:end-1,2:end-1))/dx^2+(Nutr_temp(2:end-1,1:end-2)+Nutr_temp(2:end-1,3:end)-2*Nutr_temp(2:end-1,2:end-1))/dy^2)-nutr_consump_temp(2:end-1,2:end-1);
    Phag_temp_c(2:end-1,2:end-1)=phage_product_temp(2:end-1,2:end-1);
  
    
 
    Cell_den_S_temp_c(1,:)= Cell_den_S_temp_c(2,:);  Cell_den_S_temp_c(:,1)= Cell_den_S_temp_c(:,2); 
    Cell_den_S_temp_c(end,:)= Cell_den_S_temp_c(end-1,:);   Cell_den_S_temp_c(:,end)= Cell_den_S_temp_c(:,end-1); 
    Cell_den_I_temp_c(1,:)= Cell_den_I_temp_c(2,:);  Cell_den_I_temp_c(:,1)= Cell_den_I_temp_c(:,2); 
    Cell_den_I_temp_c(end,:)= Cell_den_I_temp_c(end-1,:);   Cell_den_I_temp_c(:,end)= Cell_den_I_temp_c(:,end-1); 
    Cell_den_R_temp_c(1,:)= Cell_den_R_temp_c(2,:);  Cell_den_R_temp_c(:,1)= Cell_den_R_temp_c(:,2); 
    Cell_den_R_temp_c(end,:)= Cell_den_R_temp_c(end-1,:);   Cell_den_R_temp_c(:,end)= Cell_den_R_temp_c(:,end-1); 
    
    Attr_temp_c(1,:)=Attr_temp_c(2,:); Attr_temp_c(:,1)=Attr_temp_c(:,2); 
    Attr_temp_c(end,:)=Attr_temp_c(end-1,:);  Attr_temp_c(:,end)=Attr_temp_c(:,end-1); 
    Nutr_temp_c(1,:)=Nutr_temp_c(2,:); Nutr_temp_c(:,1)=Nutr_temp_c(:,2); 
    Nutr_temp_c(end,:)=Nutr_temp_c(end-1,:);  Nutr_temp_c(:,end)=Nutr_temp_c(:,end-1); 
    Phag_temp_c(1,:)=Phag_temp_c(2,:); Phag_temp_c(:,1)=Phag_temp_c(:,2); 
    Phag_temp_c(end,:)=Phag_temp_c(end-1,:);  Phag_temp_c(:,end)=Phag_temp_c(:,end-1); 
    
    %%
    Cell_den_S_temp=Cell_den_S+Cell_den_S_temp_c*dt;
    Cell_den_S_temp(Cell_den_S_temp<=1E-10)=0;
    Cell_den_I_temp=Cell_den_I+Cell_den_I_temp_c*dt;
    Cell_den_I_temp(Cell_den_I_temp<=1E-10)=0;
    Cell_den_R_temp=Cell_den_R+Cell_den_R_temp_c*dt;
    Cell_den_R_temp(Cell_den_R_temp<=1E-10)=0;
    Attr_temp=Attr+Attr_temp_c*dt;
    Attr_temp(Attr_temp<=0)=0;
    Nutr_temp=Nutr+Nutr_temp_c*dt; 
    Nutr_temp(Nutr_temp<=0)=0;
    Phag_temp=Phag+Phag_temp_c*dt; 
    Phag_temp(Phag_temp<=1E-10)=0;
    
    
    lambda_np_temp=lamnp0./(1+Nk_2./Nutr_temp./Nutr_temp);
    kappa_temp=kap0./(1+PT0./Phag_temp./Phag_temp);
    Gamma_temp=G0./(1+Sk./Attr_temp);
    
    S_growth_temp=lambda_np_temp.*Cell_den_S_temp-kappa_temp.*lambda_np_temp.*Cell_den_S_temp;
    I_growth_temp=beta1*lambda_np_temp.*Cell_den_I_temp+kappa_temp.*lambda_np_temp.*Cell_den_S_temp-theta*Cell_den_I_temp;
    R_growth_temp=beta2*lambda_np_temp.*Cell_den_R_temp+theta*Cell_den_I_temp;
    
    nutr_consump_temp=lambda_np_temp.*(Cell_den_S_temp+Cell_den_I_temp+Cell_den_R_temp)./Yn;
    attr_consump_temp=Gamma_temp.*(Cell_den_S_temp+Cell_den_I_temp+Cell_den_R_temp);
    phage_product_temp=eta*((1-beta1)*lambda_np_temp.*Cell_den_I_temp+(1-beta2)*lambda_np_temp.*Cell_den_R_temp);
    
    
    f=chi.*log((1+ Attr_temp/K1)./(1+ Attr_temp/K2));% free energy
    g1_S(2:end-1,:)=Cell_den_S_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_S(:,2:end-1)=Cell_den_S_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_S(1,:)=-g1_S(2,:);g1_S(end,:)=-g1_S(end-1,:);
    g2_S(:,1)=-g2_S(:,2);g2_S(:,end)=-g2_S(:,end-1); 
    
    g1_I(2:end-1,:)=Cell_den_I_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_I(:,2:end-1)=Cell_den_I_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_I(1,:)=-g1_I(2,:);g1_I(end,:)=-g1_I(end-1,:);
    g2_I(:,1)=-g2_I(:,2);g2_I(:,end)=-g2_I(:,end-1); 
    
    g1_R(2:end-1,:)=Cell_den_R_temp(2:end-1,:).*(f(3:end,:)-f(1:end-2,:))/2/dx;
    g2_R(:,2:end-1)=Cell_den_R_temp(:,2:end-1).*(f(:,3:end)-f(:,1:end-2))/2/dy;
    g1_R(1,:)=-g1_R(2,:);g1_R(end,:)=-g1_R(end-1,:);
    g2_R(:,1)=-g2_R(:,2);g2_R(:,end)=-g2_R(:,end-1); 


    
    Cell_den_S_temp_d(2:end-1,2:end-1)=Diff.*((Cell_den_S_temp(1:end-2,2:end-1)+Cell_den_S_temp(3:end,2:end-1)-2*Cell_den_S_temp(2:end-1,2:end-1))/dx^2+(Cell_den_S_temp(2:end-1,1:end-2)+Cell_den_S_temp(2:end-1,3:end)-2*Cell_den_S_temp(2:end-1,2:end-1))/dy^2)-(g1_S(3:end,2:end-1)-g1_S(1:end-2,2:end-1))/2/dx-(g2_S(2:end-1,3:end)-g2_S(2:end-1,1:end-2))/2/dy+S_growth_temp(2:end-1,2:end-1);
    Cell_den_I_temp_d(2:end-1,2:end-1)=Diff.*((Cell_den_I_temp(1:end-2,2:end-1)+Cell_den_I_temp(3:end,2:end-1)-2*Cell_den_I_temp(2:end-1,2:end-1))/dx^2+(Cell_den_I_temp(2:end-1,1:end-2)+Cell_den_I_temp(2:end-1,3:end)-2*Cell_den_I_temp(2:end-1,2:end-1))/dy^2)-(g1_I(3:end,2:end-1)-g1_I(1:end-2,2:end-1))/2/dx-(g2_I(2:end-1,3:end)-g2_I(2:end-1,1:end-2))/2/dy+I_growth_temp(2:end-1,2:end-1);
    Cell_den_R_temp_d(2:end-1,2:end-1)=Diff.*((Cell_den_R_temp(1:end-2,2:end-1)+Cell_den_R_temp(3:end,2:end-1)-2*Cell_den_R_temp(2:end-1,2:end-1))/dx^2+(Cell_den_R_temp(2:end-1,1:end-2)+Cell_den_R_temp(2:end-1,3:end)-2*Cell_den_R_temp(2:end-1,2:end-1))/dy^2)-(g1_R(3:end,2:end-1)-g1_R(1:end-2,2:end-1))/2/dx-(g2_R(2:end-1,3:end)-g2_R(2:end-1,1:end-2))/2/dy+R_growth_temp(2:end-1,2:end-1);
    
    Attr_temp_d(2:end-1,2:end-1)=Da*((Attr_temp(1:end-2,2:end-1)+Attr_temp(3:end,2:end-1)-2*Attr_temp(2:end-1,2:end-1))/dx^2+(Attr_temp(2:end-1,1:end-2)+Attr_temp(2:end-1,3:end)-2*Attr_temp(2:end-1,2:end-1))/dy^2)-attr_consump_temp(2:end-1,2:end-1);
    Nutr_temp_d(2:end-1,2:end-1)=Dn*((Nutr_temp(1:end-2,2:end-1)+Nutr_temp(3:end,2:end-1)-2*Nutr_temp(2:end-1,2:end-1))/dx^2+(Nutr_temp(2:end-1,1:end-2)+Nutr_temp(2:end-1,3:end)-2*Nutr_temp(2:end-1,2:end-1))/dy^2)-nutr_consump_temp(2:end-1,2:end-1);
    Phag_temp_d(2:end-1,2:end-1)=phage_product_temp(2:end-1,2:end-1);
    
 
    Cell_den_S_temp_d(1,:)= Cell_den_S_temp_d(2,:);  Cell_den_S_temp_d(:,1)= Cell_den_S_temp_d(:,2); 
    Cell_den_S_temp_d(end,:)= Cell_den_S_temp_d(end-1,:);   Cell_den_S_temp_d(:,end)= Cell_den_S_temp_d(:,end-1); 
    Cell_den_I_temp_d(1,:)= Cell_den_I_temp_d(2,:);  Cell_den_I_temp_d(:,1)= Cell_den_I_temp_d(:,2); 
    Cell_den_I_temp_d(end,:)= Cell_den_I_temp_d(end-1,:);   Cell_den_I_temp_d(:,end)= Cell_den_I_temp_d(:,end-1); 
    Cell_den_R_temp_d(1,:)= Cell_den_R_temp_d(2,:);  Cell_den_R_temp_d(:,1)= Cell_den_R_temp_d(:,2); 
    Cell_den_R_temp_d(end,:)= Cell_den_R_temp_d(end-1,:);   Cell_den_R_temp_d(:,end)= Cell_den_R_temp_d(:,end-1); 
       
    Attr_temp_d(1,:)=Attr_temp_d(2,:); Attr_temp_d(:,1)=Attr_temp_d(:,2); 
    Attr_temp_d(end,:)=Attr_temp_d(end-1,:);  Attr_temp_d(:,end)=Attr_temp_d(:,end-1); 
    Nutr_temp_d(1,:)=Nutr_temp_d(2,:); Nutr_temp_d(:,1)=Nutr_temp_d(:,2); 
    Nutr_temp_d(end,:)=Nutr_temp_d(end-1,:);  Nutr_temp_d(:,end)=Nutr_temp_d(:,end-1); 
    Phag_temp_d(1,:)=Phag_temp_d(2,:); Phag_temp_d(:,1)=Phag_temp_d(:,2); 
    Phag_temp_d(end,:)=Phag_temp_d(end-1,:);  Phag_temp_d(:,end)=Phag_temp_d(:,end-1);     
    %%
    Cell_den_S=Cell_den_S+1/6*(Cell_den_S_temp_a+2*Cell_den_S_temp_b+2*Cell_den_S_temp_c+Cell_den_S_temp_d)*dt;
    Cell_den_I=Cell_den_I+1/6*(Cell_den_I_temp_a+2*Cell_den_I_temp_b+2*Cell_den_I_temp_c+Cell_den_I_temp_d)*dt;    
    Cell_den_R=Cell_den_R+1/6*(Cell_den_R_temp_a+2*Cell_den_R_temp_b+2*Cell_den_R_temp_c+Cell_den_R_temp_d)*dt;  
    Attr=Attr+1/6*(Attr_temp_a+2*Attr_temp_b+2*Attr_temp_c+Attr_temp_d)*dt;
    Nutr=Nutr+1/6*(Nutr_temp_a+2*Nutr_temp_b+2*Nutr_temp_c+Nutr_temp_d)*dt;
    Phag=Phag+1/6*(Phag_temp_a+2*Phag_temp_b+2*Phag_temp_c+Phag_temp_d)*dt;
    
   Cell_den_S(1,:)=Cell_den_S(2,:); Cell_den_S(:,1)=Cell_den_S(:,2); 
   Cell_den_S(end,:)=Cell_den_S(end-1,:);  Cell_den_S(:,end)=Cell_den_S(:,end-1);   
   Cell_den_S(Cell_den_S<=1E-10)=0; 
   
   Cell_den_I(1,:)=Cell_den_I(2,:); Cell_den_I(:,1)=Cell_den_I(:,2); 
   Cell_den_I(end,:)=Cell_den_I(end-1,:);  Cell_den_I(:,end)=Cell_den_I(:,end-1);   
   Cell_den_I(Cell_den_I<=1E-10)=0;  
   
   Cell_den_R(1,:)=Cell_den_R(2,:); Cell_den_R(:,1)=Cell_den_R(:,2); 
   Cell_den_R(end,:)=Cell_den_R(end-1,:);  Cell_den_R(:,end)=Cell_den_R(:,end-1);   
   Cell_den_R(Cell_den_R<=1E-10)=0; 
    
    
    Attr(1,:)=Attr(2,:); Attr(:,1)=Attr(:,2); 
    Attr(end,:)=Attr(end-1,:);  Attr(:,end)=Attr(:,end-1);
    Attr(Attr<=0)=0; 
    Nutr(1,:)=Nutr(2,:); Nutr(:,1)=Nutr(:,2); 
    Nutr(end,:)=Nutr(end-1,:);  Nutr(:,end)=Nutr(:,end-1); 
    Nutr(Nutr<=0)=0;
    Phag(1,:)=Phag(2,:); Phag(:,1)=Phag(:,2); 
    Phag(end,:)=Phag(end-1,:);  Phag(:,end)=Phag(:,end-1);  
    Phag(Phag<=1E-10)=0; 
    
    

       
    if t/3600==round(t/3600)
    save(strrep(name,'.mat',['_' num2str(t/3600) '.mat']));
    end   
    
end

toc

end
