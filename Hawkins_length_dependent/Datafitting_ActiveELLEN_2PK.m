function f=Datafitting_ActiveELLEN_2PK(para)
global active cmax
%%
%parapassive=[3.2790    0.0010    0.1418    4.6476];
%i=0; t0=0.01; tr=0.5;
%lmax=1.4;inte=0;lamd=lmax;pmax=104.72;k=0.25;vmax=193;kr=0.005;
k=0;

 beta=3;
 fc=0.5+1/pi*atan(beta*log(cmax));
 
 for i=1:size(active,1)
     k=k+1;
     lamd=active(i,1);

     %lamax=(para(1)*lamd+para(2));
     %lamax=para(5);
     lamax=(para(1)*lamd^2+para(2)*lamd+para(3));

     ps=(0-1)/(0-lamax);

     ta=ps*lamax/(1+fc*(ps-1));

     dd=(1-ta)*(9.25*lamd-1);
     lamdf=lamd*9.25-ta*(9.25*lamd-1);

     po=para(4)*(lamdf-1)*lamdf^2;
     
     f(k)=po-active(i,2);

 end
 
end
