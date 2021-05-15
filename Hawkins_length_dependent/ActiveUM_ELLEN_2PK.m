clear all;close all;

 global active cmax

 cmax=1;

 
 %dddd
for ii=1:20
load('EXPdataCau');

para0=1*rand(1,5);
%para0=[ 0.5281   -0.5584    1.4179  387.8688    0.4167];

%boundary values
Lb=[-10,-10,-10,10,0.6];
Ub=[10,10,10,100,1];

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_ActiveELLEN_2PK(para),para0,Lb,Ub,options);
para0=para;

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_ActiveELLEN_2PK(para),para0,Lb,Ub,options);
para;
abadd(ii)=resnorm;
dadd(ii,:)=para;

end

eges=find(abadd==min(abadd));
para=real(dadd(eges,:))
 
 %
  beta=3;
 fc=0.5+1/pi*atan(beta*log(cmax));
 
 for i=1:size(active,1)
     lamd=active(i,1);
 
     %lamax=(para(1)*lamd+para(2));
     %lamax=para(5);
     lamax=(para(1)*lamd^2+para(2)*lamd+para(3));
     
     ps=(0-1)/(0-lamax);
     
     ta=ps*lamax/(1+fc*(ps-1));

    
     dd=(1-ta)*(9.25*lamd-1);

     lamdf=lamd*9.25-ta*(9.25*lamd-1);
     po(i,2)=para(4)*(lamdf-1)*lamdf^2;
     
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
 k=0;parap=[3.2790    0.0010    0.1418    4.6476];
for i=1:size(total,1)
    lamd=total(i,1);
    
     %lamax=(para(1)*lamd+para(2));
     %lamax=para(5);
     lamax=(para(1)*lamd^2+para(2)*lamd+para(3));
     
     ps=(0-1)/(0-lamax);
     %ta=lamax/((1-fc)*lamax+fc);
     ta=ps*lamax/(1+fc*(ps-1));
     dd=(1-ta)*(9.25*lamd-1);
     %dd=lamd*(1.54-ta*1.375)-1;
     %dd=(1-ta)*(2.85-2.2*lamd);
     %dd=(1-ta)*1.65;
     
    po(i,1)=lamd;
    
    lamdf=lamd*9.25-ta*(9.25*lamd-1);
     po(i,2)=para(4)*(lamdf-1)*lamdf^2;
    
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
 resom=sum(res)/size(total,1)


figure(3)
hold on

plot(total(:,1),total(:,2))
plot(active(:,1),active(:,2))
plot(passive(:,1),passive(:,2))
plot(pt(:,1),pt(:,2))
plot(po(:,1),po(:,2))
plot(pa(:,1),pa(:,2))
    
