function f=Datafitting(para)
global dataT
tspan = [0 1000];    
y1_0 = 1;    
y2_0 = 0;    
y3_0 = 0;    
[T,Y] = ode15s(@oscfit,tspan,[y1_0 y2_0 y3_0 para(1:8)]);   


xiT=0:10:1000;
yispl=interp1(T,Y(:,3),xiT,'spline');

 beta=3;

 r1=0.1475; r2=-0.3500; r3=1.1738;
 rat=para(9);
 r1=r1*rat; r2=r2*rat; r3=r3*rat;
 
 eta=9.25;k1=para(10);
 lamd_a=2/1.65;lamd1=lamd_a/1*(1-1/eta)+1/eta;
 lamdac=1;

 for i=1:length(xiT)
     
     %lamd=lamd_a/lamdac*(1-1/eta)+1/eta;
     lamd=lamd1;
     cbar=real(Y(i,3));
     fc=1+2/pi*atan(beta*log(cbar));
     lmax=(r1*lamd^2+r2*lamd+r3);
     xi=1/lmax;
     lamdac=real(xi*lmax/(1+fc*(xi-1)));
     
     po=k1*(1-lamdac)*(eta*lamd-1);
     
     
     f(i)=po-dataT(i);

 end


return
end