
function []=phage_SIR_liquid_function(name,phage_init)
%% define parameters
dt=10;                  %time step in unit of seconds
Time=60*60*10;         %simulated time;
T=round(Time/dt);       %time steps

Nk=0.4;                    % nutrient consumption half threshold mM


Yn=0.064;                %0.064OD/mM
Nutr_init=28;            % initial nutrient concentration 40mM
bacter_init=0.005;
% phage_init=0.01;

lamnp0=log(2)/20/60;
kap0=0.7;

P0=0.004;

alpha=700;
beta1=0.1;
beta2=0.9;
theta=1E-4;
xi=5E-4;
%% initial conditions
Cell_den_1=bacter_init;
Cell_den_2=0;
Cell_den_3=0;
Nutr=Nutr_init;
Phag=phage_init;


tic

%%
for t=1:T

    
    Cell_den_1_temp=Cell_den_1;
    Cell_den_1_temp(Cell_den_1_temp<=0)=0;
    Cell_den_2_temp=Cell_den_2;
    Cell_den_2_temp(Cell_den_2_temp<=0)=0;
    Cell_den_3_temp=Cell_den_3;
    Cell_den_3_temp(Cell_den_3_temp<=0)=0;
    Nutr_temp=Nutr; 
    Nutr_temp(Nutr_temp<=0)=0;
    Phag_temp=Phag; 
    Phag_temp(Phag_temp<=0)=0;
    
    lambda_np_temp=lamnp0./(1+Nk/Nutr_temp);
    
    kappa_temp=kap0;
    
    Cell_den_1_temp_a=lambda_np_temp*Cell_den_1_temp-kappa_temp.*Phag_temp.*lambda_np_temp*Cell_den_1_temp;
    Cell_den_2_temp_a=beta1*lambda_np_temp*Cell_den_2_temp+kappa_temp.*Phag_temp.*lambda_np_temp*Cell_den_1_temp-theta*Cell_den_2_temp;
    Cell_den_3_temp_a=beta2*lambda_np_temp*Cell_den_3_temp+theta*Cell_den_2_temp;
    Nutr_temp_a=-lambda_np_temp.*(Cell_den_1_temp+Cell_den_2_temp+Cell_den_3_temp)./Yn;
    Phag_temp_a=alpha*((1-beta1)*lambda_np_temp*Cell_den_2_temp+(1-beta2)*lambda_np_temp*Cell_den_3_temp)-xi*Cell_den_1_temp;

    %%
    Cell_den_1_temp=Cell_den_1+Cell_den_1_temp_a*dt/2;
    Cell_den_1_temp(Cell_den_1_temp<=0)=0;
    Cell_den_2_temp=Cell_den_2+Cell_den_2_temp_a*dt/2;
    Cell_den_2_temp(Cell_den_2_temp<=0)=0;
    Cell_den_3_temp=Cell_den_3+Cell_den_3_temp_a*dt/2;
    Cell_den_3_temp(Cell_den_3_temp<=0)=0;
    Nutr_temp=Nutr+Nutr_temp_a*dt/2; 
    Nutr_temp(Nutr_temp<=0)=0;
    Phag_temp=Phag+Phag_temp_a*dt/2; 
    Phag_temp(Phag_temp<=0)=0;
    
    
    lambda_np_temp=lamnp0./(1+Nk/Nutr_temp);
    
    kappa_temp=kap0;  
    
    Cell_den_1_temp_b=lambda_np_temp*Cell_den_1_temp-kappa_temp.*Phag_temp.*lambda_np_temp*Cell_den_1_temp;
    Cell_den_2_temp_b=beta1*lambda_np_temp*Cell_den_2_temp+kappa_temp.*Phag_temp.*lambda_np_temp*Cell_den_1_temp-theta*Cell_den_2_temp;
    Cell_den_3_temp_b=beta2*lambda_np_temp*Cell_den_3_temp+theta*Cell_den_2_temp;
    Nutr_temp_b=-lambda_np_temp.*(Cell_den_1_temp+Cell_den_2_temp+Cell_den_3_temp)./Yn;
    Phag_temp_b=alpha*((1-beta1)*lambda_np_temp*Cell_den_2_temp+(1-beta2)*lambda_np_temp*Cell_den_3_temp)-xi*Cell_den_1_temp;
    

    %%
    Cell_den_1_temp=Cell_den_1+Cell_den_1_temp_b*dt/2;
    Cell_den_1_temp(Cell_den_1_temp<=0)=0;
    Cell_den_2_temp=Cell_den_2+Cell_den_2_temp_b*dt/2;
    Cell_den_2_temp(Cell_den_2_temp<=0)=0;
    Cell_den_3_temp=Cell_den_3+Cell_den_3_temp_b*dt/2;
    Cell_den_3_temp(Cell_den_3_temp<=0)=0;
    
    Nutr_temp=Nutr+Nutr_temp_b*dt/2; 
    Nutr_temp(Nutr_temp<=0)=0;
    Phag_temp=Phag+Phag_temp_b*dt/2; 
    Phag_temp(Phag_temp<=0)=0;
    
    
    lambda_np_temp=lamnp0./(1+Nk/Nutr_temp);
    
    kappa_temp=kap0; 
    
    
    Cell_den_1_temp_c=lambda_np_temp*Cell_den_1_temp-kappa_temp.*Phag_temp.*lambda_np_temp*Cell_den_1_temp;
    Cell_den_2_temp_c=beta1*lambda_np_temp*Cell_den_2_temp+kappa_temp.*Phag_temp.*lambda_np_temp*Cell_den_1_temp-theta*Cell_den_2_temp;
    Cell_den_3_temp_c=beta2*lambda_np_temp*Cell_den_3_temp+theta*Cell_den_2_temp;
    Nutr_temp_c=-lambda_np_temp.*(Cell_den_1_temp+Cell_den_2_temp+Cell_den_3_temp)./Yn;
    Phag_temp_c=alpha*((1-beta1)*lambda_np_temp*Cell_den_2_temp+(1-beta2)*lambda_np_temp*Cell_den_3_temp)-xi*Cell_den_1_temp; 
    %%
    Cell_den_1_temp=Cell_den_1+Cell_den_1_temp_c*dt;
    Cell_den_1_temp(Cell_den_1_temp<=0)=0;
    Cell_den_2_temp=Cell_den_2+Cell_den_2_temp_c*dt;
    Cell_den_2_temp(Cell_den_2_temp<=0)=0;
    Cell_den_3_temp=Cell_den_3+Cell_den_3_temp_c*dt;
    Cell_den_3_temp(Cell_den_3_temp<=0)=0;
    Nutr_temp=Nutr+Nutr_temp_c*dt; 
    Nutr_temp(Nutr_temp<=0)=0;
    Phag_temp=Phag+Phag_temp_c*dt; 
    Phag_temp(Phag_temp<=0)=0;
    
    
    lambda_np_temp=lamnp0./(1+Nk/Nutr_temp);
    
    kappa_temp=kap0; 
    
    Cell_den_1_temp_d=lambda_np_temp*Cell_den_1_temp-kappa_temp.*Phag_temp.*lambda_np_temp*Cell_den_1_temp;
    Cell_den_2_temp_d=beta1*lambda_np_temp*Cell_den_2_temp+kappa_temp.*Phag_temp.*lambda_np_temp*Cell_den_1_temp-theta*Cell_den_2_temp;
    Cell_den_3_temp_d=beta2*lambda_np_temp*Cell_den_3_temp+theta*Cell_den_2_temp;
    Nutr_temp_d=-lambda_np_temp.*(Cell_den_1_temp+Cell_den_2_temp+Cell_den_3_temp)./Yn;
    Phag_temp_d=alpha*((1-beta1)*lambda_np_temp*Cell_den_2_temp+(1-beta2)*lambda_np_temp*Cell_den_3_temp)-xi*Cell_den_1_temp;
    
    %%
    Cell_den_1=Cell_den_1+1/6*(Cell_den_1_temp_a+2*Cell_den_1_temp_b+2*Cell_den_1_temp_c+Cell_den_1_temp_d)*dt;
    Cell_den_2=Cell_den_2+1/6*(Cell_den_2_temp_a+2*Cell_den_2_temp_b+2*Cell_den_2_temp_c+Cell_den_2_temp_d)*dt;
    Cell_den_3=Cell_den_3+1/6*(Cell_den_3_temp_a+2*Cell_den_3_temp_b+2*Cell_den_3_temp_c+Cell_den_3_temp_d)*dt;
    
    Nutr=Nutr+1/6*(Nutr_temp_a+2*Nutr_temp_b+2*Nutr_temp_c+Nutr_temp_d)*dt;
    Phag=Phag+1/6*(Phag_temp_a+2*Phag_temp_b+2*Phag_temp_c+Phag_temp_d)*dt;
    
    Cell_den_1(Cell_den_1<=0)=0; 
    Cell_den_2(Cell_den_2<=0)=0; 
    Cell_den_3(Cell_den_3<=0)=0; 
    Nutr(Nutr<=0)=0;
    Phag(Phag<=0)=0; 
    
    Cell_1_Save(t)=Cell_den_1;
    Cell_2_Save(t)=Cell_den_2;
    Cell_3_Save(t)=Cell_den_3;
    Cell_Save(t)=Cell_den_1+Cell_den_2+Cell_den_3;
    Nutr_Save(t)=Nutr;
    Phag_Save(t)=Phag; 
    
    
    if t/360==round(t/360)
    save(strrep(name,'.mat',['_' num2str(t/360) '.mat']));
    end  
 

end


toc
end