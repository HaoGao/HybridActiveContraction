function f=Datafitting_ActiveELLEN(para)
global active cmax
%%
%parapassive=[3.2790    0.0010    0.1418    4.6476];
%i=0; t0=0.01; tr=0.5;
%lmax=1.4;inte=0;lamd=lmax;pmax=104.72;k=0.25;vmax=193;kr=0.005;
k=0;

 beta=3; para(5)=9.25;
 %fc=(0.5+1/pi*atan(beta*log(cmax)));
 fc=1+2/pi*atan(beta*log(cmax));
 
 for i=1:size(active,1)
     k=k+1;
     lamd=active(i,1);
     %lamax=0.5;
     %if lamd>1
     %lamax=(1+para(2)*(lamd-1));
     %else
     %    lamax=(1+para(2)*(1-1));
     %end
     %lamax=sin(-para(1)*lamd+para(2));
     %para(3)=0.89+para(2)*para(2)/4/para(1);
      lamax=(para(1)*lamd^2+para(2)*lamd+para(3));
     %if lamax>0 & lamax <1
     ps=(0-1)/(0-lamax);
     %ta=lamax/((1-fc)*lamax+fc);
     ta=ps*lamax/(1+fc*(ps-1));
    % ta=lamax;
     
     %dd=(1-para(3))*ta+(para(3)*lamd-1);
     dd=(1-ta)*(para(5)*lamd-1);
     %dd=(1-ta)*1.65;
     %dd=lamd*(1.54-ta*1.375)-1;
     %ps=(0-1)/(0-lamax);
     %ta=ps*lamax/(1+fc*(ps-1));
     %dd=(1-ta)*(3.85-2.2*lamd);
     
     %dd=(1-ta)*(para(3)*lamd-1)*(-para(1)*lamd+para(2));;

     %po=para(4)*dd;
     po=250*dd;
     %if lamax>1
     %    f(k)=(po-active(i,2))*10000;
     %else
     f(k)=po-active(i,2);
     %end
     %else
     %    f(k)=10000000000000;
     %end
 end
 
end
