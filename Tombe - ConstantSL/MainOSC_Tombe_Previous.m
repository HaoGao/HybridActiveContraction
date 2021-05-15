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
 %r1=1.2749;r2=-2.8657;r3=2.3467;
 r1=0.5017; r2=-1.1265; r3=1.5332;
 eta=9.25;k1=130;
 lamd_a=2/1.65;lamd1=lamd_a/1*(1-1/eta)+1/eta;
 lamdac=1;
 
 parap=[3.2790    0.0010    0.1418    4.6476];
 
 for i=1:length(T)
     
     %lamd=lamd_a/lamdac*(1-1/eta)+1/eta;
     lamd=lamd1;
     
     
     cbar=real(Y(i,3));
     fc=0.5+1/pi*atan(beta*log(cbar));
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
 plot(T,lamdSLE)
 plot(T,lamdSLM)
 plot(T,lamd_la)
 plot(T,lamdSLf)
 
 
 %load('Tombe_ML');
 %plot(Tombe_ML(:,1),Tombe_ML(:,2)/Tombe_ML(1,2))

 figure(7)
 hold on
 plot(T,po)
% plot(T,pp-pp(1))
% plot(T,pt)
 
 
 load('Tombe_MLForce');
 plot(Tombe_MLForce(:,1),Tombe_MLForce(:,2))
 
% load('GU-Force');
% plot(force(:,1)-50,force(:,2)*50)

     