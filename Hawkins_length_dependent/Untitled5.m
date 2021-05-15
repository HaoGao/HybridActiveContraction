clear all;%close all;
%
 global active cmax
 
 %
 %tspan = [0 500];   Environmental Physiology of Animals
 %y1_0 = 1;
 %y2_0 = 0;
 %y3_0 = 0;
 %[T,Y] = ode15s(@osc,tspan,[y1_0 y2_0 y3_0]);
 %figure(1)
 %hold on
 %plot(T,Y(:,1),'-')
 %plot(T,Y(:,2),'-')
 %plot(T,Y(:,3),'-')
 
 %for i=1:length(T)
 %    lads(i)=1/(1+Y(i,3)*0.9);
 %end
 %plot(T,lads)
 cmax=1;

 
 %dddd

load('EXPdataCau');
%ii=16;
ks1=16*10;

para0=1*rand(1,5);
%para0=[ 0.5281   -0.5584    1.4179  387.8688    0.4167];

%boundary values
Lb=[-100,-100,-100,ks1,0];
Ub=[100,100,100,500,1000];

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_ActiveELLEN(para),para0,Lb,Ub,options);
para0=para;

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_ActiveELLEN(para),para0,Lb,Ub,options);
para






 
para(5)=9.25;
 %
  beta=3;
 fc=(0.5+1/pi*atan(beta*log(cmax)))*2;
 
 for i=1:size(active,1)
     lamd=active(i,1);
 %    if lamd>1
 %    lamax=(1+para(2)*(lamd-1));
 %    else
 %        lamax=(1+para(2)*(1-1));
 %    end
     lamax=(para(1)*lamd^2+para(2)*lamd+para(3));
     ps=(0-1)/(0-lamax);
     %ta=lamax/((1-fc)*lamax+fc);
     ta=ps*lamax/(1+fc*(ps-1));

     
     %dd=(1-para(3))*ta+(para(3)*lamd-1); 
     dd=(1-ta)*(para(5)*lamd-1);
     %dd=lamd*(1.54-ta*1.375)-1;
     %dd=(1-ta)*(3.85-2.2*lamd);
     % dd=(1-ta)*1.65;

     po(i,2)=para(4)*dd;
     
     po(i,1)=lamd;
     tl(i,1)=lamd;tl(i,2)=ta;tl(i,3)=lamax;tl(i,4)=dd;
 end
 
 figure(2)
 plot(tl(:,1),tl(:,2))
 hold on
 plot(tl(:,1),tl(:,3))
 plot(tl(:,1),tl(:,4))
 
 min(tl(:,3))
  max(tl(:,3))
 %
 k=0;%parap=[3.2790    0.0010    0.1418    4.6476];
 parap=[3.6039    0.1000    0.2063    4.6250 ];
for i=1:size(total,1)
    lamd=total(i,1);
    
     lamax=(para(1)*lamd^2+para(2)*lamd+para(3));
     ps=(0-1)/(0-lamax);
     %ta=lamax/((1-fc)*lamax+fc);
     ta=ps*lamax/(1+fc*(ps-1));
     dd=(1-ta)*(para(5)*lamd-1);
     %dd=lamd*(1.54-ta*1.375)-1;
     %dd=(1-ta)*(2.85-2.2*lamd);
     %dd=(1-ta)*1.65;
     
    po(i,1)=lamd;
     po(i,2)=para(4)*dd;
    
    if lamd>1
        k=k+1;
        [sigma_fs]=InteractionUM(lamd,parap);
        pa(k,1)=lamd;
        pa(k,2)=sigma_fs;
        pt(i,2)=po(i,2)+sigma_fs;
    else
         pt(i,2)=po(i,2);
    end
    res(i)=(pt(i,2)-total(i,2))^2;
    pt(i,1)=lamd;
    
end



figure(3)
hold on

plot(total(:,1),total(:,2))
plot(active(:,1),active(:,2))
plot(passive(:,1),passive(:,2))
plot(pt(:,1),pt(:,2))
plot(po(:,1),po(:,2))
plot(pa(:,1),pa(:,2))
    




