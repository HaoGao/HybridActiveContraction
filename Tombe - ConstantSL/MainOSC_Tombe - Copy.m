clear all; 
%close all;

tspan = [0 1000];    
y1_0 = 1;    
y2_0 = 0;    
y3_0 = 0;    
[T,Y] = ode15s(@osc,tspan,[y1_0 y2_0 y3_0]);   

%figure(1)
% plot(T,Y(:,1),'-') 
% hold on
 
% plot(T,Y(:,2),'-') 
% plot(T,Y(:,3),'-') 
 
 beta=3;
 r1=1.2749;r2=-2.8657;r3=2.3467;
 eta=9.25;k1=45;
 lamd_a=2.2/1.85;lamd1=lamd_a/1*(1-1/eta)+1/eta;
 lamdac=1;
 for i=1:length(T)
     
     lamd=lamd_a/lamdac*(1-1/eta)+1/eta;
     lamdML(i)=lamd/lamd1;
     
     cbar=real(Y(i,3));
     fc=0.5+1/pi*atan(beta*log(cbar));
     lmax=(r1*lamd^2+r2*lamd+r3)*0.96;
     xi=1/lmax;
     lamdac=real(xi*lmax/(1+fc*(xi-1)));
     
     po(i)=k1*(1-lamdac)*(eta*lamd-1);

 end
 figure(6)
 hold on
 %plot(T,lamd_a)
 plot(T,lamdML)
 
 figure(7)
 hold on
 plot(T,po)
 %plot(T,lamdML)
 
% load('GU-Force');
% plot(force(:,1)-50,force(:,2)*50)

     