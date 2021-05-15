clear all; 
global dataT
%close all;
load('Tombe_MLForce');

xiT2=0:10:1000;
dataT=interp1(Tombe_MLForce(:,1),Tombe_MLForce(:,2),xiT2,'spline');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

para0=[0.01 0.15 1 0.04 0.3 0.002  0.015 0.015  1.007 250];

%boundary values
Lb=[0.0080    0.1200    0.8000    0.0320    0.2400    0.0016    0.0120    0.0120    0.9  200];
Ub=[0.0150    0.2250    1.5000    0.0600    0.4500    0.0030    0.0225    0.0225    1.1  375];

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting(para),para0,Lb,Ub,options);
para0=para;
[para,resnorm,residual]=lsqnonlin(@(para)Datafitting(para),para0,Lb,Ub,options);
para0=para;
[para,resnorm,residual]=lsqnonlin(@(para)Datafitting(para),para0,Lb,Ub,options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tspan = [0 1000];    
y1_0 = 1;    
y2_0 = 0;    
y3_0 = 0;    
[T,Y] = ode15s(@oscfit,tspan,[y1_0 y2_0 y3_0 para(1:8)]);   

figure(1)
 plot(T/1000,Y(:,1),'-') 
 hold on
 
 plot(T/1000,Y(:,2),'-') 
 plot(T/1000,Y(:,3),'-') 
 
 beta=3;

 r1=0.1475; r2=-0.3500; r3=1.1738;
 rat=para(9);
 r1=r1*rat; r2=r2*rat; r3=r3*rat;
 
 eta=9.25;k1=para(10);
 lamd_a=2/1.65;lamd1=lamd_a/1*(1-1/eta)+1/eta;
 lamdac=1;
 
 for i=1:length(T)
     
     %lamd=lamd_a/lamdac*(1-1/eta)+1/eta;
     lamd=lamd1;
     
     
     cbar=real(Y(i,3));
     fc=1+2/pi*atan(beta*log(cbar));
     lmax=(r1*lamd^2+r2*lamd+r3);
     xi=1/lmax;
     lamdac=real(xi*lmax/(1+fc*(xi-1)));
     lamdSLE(i)=lamdac;
     lamdSLM(i)=lamd_a;
     lamd_la(i)=lamdac*lamd_a;
     lamdSLf(i)=(lamd1*1.85-lamd_la(i)*1.65)/0.2;
     
     po(i)=k1*(1-lamdac)*(eta*lamd-1);
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %    [sigma_fs]=InteractionUM(lamd,parap);
    %    pp(i)=sigma_fs;
    %    pt(i)=po(i)+pp(i)-pp(1);
     
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 end
 

 figure(6)
 hold on
 %plot(T,lamd_a)
 plot(T/1000,lamdSLE)
 plot(T/1000,lamdSLM)
 plot(T/1000,lamd_la)
 plot(T/1000,lamdSLf)
 
 
 %load('Tombe_ML');
 %plot(Tombe_ML(:,1),Tombe_ML(:,2)/Tombe_ML(1,2))

 figure(7)
 hold on
 plot(T/1000,po)
% plot(T,pp-pp(1))
% plot(T,pt)
 
 
 plot(Tombe_MLForce(:,1)/1000,Tombe_MLForce(:,2))
 
% load('GU-Force');
% plot(force(:,1)-50,force(:,2)*50)

     