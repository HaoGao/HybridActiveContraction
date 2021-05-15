clear all; %close all;
global passive total active

load('EXPdataCau');
%load('lamda');

for ii=1:10

para0=1*rand(1,10)

%boundary values

Lb=[-100,-1000,1e-3,1e-3,1e-3,1e-3,   1e-3,0.5,   -50,0];
Ub=[100,5000,50,50,50,50,  50,1,   50,50];


options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_ActiveFEFA_ADD(para),para0,Lb,Ub,options);
para0=para;


[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_ActiveFEFA_ADD(para),para0,Lb,Ub,options);

abadd(ii)=resnorm;
dadd(ii,:)=para;

end

eges=find(abadd==min(abadd));
para=real(dadd(eges,:))

parap=[3.6039    0.1000    0.2063    4.6250];
    a=parap(1);b=parap(2);
    a_f=parap(3);b_f=parap(4);
    eps=para(7);
    
figure(1)
hold on

plot(total(:,1),total(:,2))
plot(active(:,1),active(:,2))
plot(passive(:,1),passive(:,2))

k=0;
 fc=1;%0.3648;
for i=1:size(total,1)
    lamd=total(i,1);
    
    F=[lamd 0 0; 0 1/sqrt(lamd) 0; 0 0 1/sqrt(lamd)];
    %lamd_a=call_osc4_plot(lamd,para);
    %lamax=(para(1)*lamd+para(2));
    %lamax=(para(8)*lamd^2+para(9)*lamd+para(10));
    lamax=para(8);
     
     ps=(0-1)/(0-lamax);
     %ta=lamax/((1-fc)*lamax+fc);
    lamd_a=ps*lamax/(1+fc*(ps-1));
    
    tads(i)=lamd_a;

%invariants
    f0=[1 0 0]';
    FEE=F+(1/lamd_a-1)*(F*f0)*f0';
    FE=FEE;%*det(FEE)^(-1/3);
    detm=det(FE);

    CE=FE'*FE;
    I_1=sum(diag(CE));
    I_4f=dot(f0,CE*f0);
    I_4_f=max(I_4f,1);
    BE=FE*FE';
    fxfE=(FE*f0)*(FE*f0)';
    
    c2=2*eps*(I_4_f-1)*fxfE;  
    sigma1=[];
    sigma1=[c2(1,1)-c2(3,3)];
    
    ac(i,2)=sum(sigma1)/detm;
    ac(i,1)=lamd;
    
    
    C=F'*F;
    I_1=sum(diag(C));
    I_4f=dot(f0,C*f0);
    I_4_f=max(I_4f,1);
    B=F*F';
    fxf=(F*f0)*(F*f0)';
%stress

    if I_4_f>1
        k=k+1;
        c1=(a*exp(b*(I_1-3)))*B;
        c2=(2*a_f*(I_4_f-1)*exp(b_f*(I_4_f-1)^2))*fxf;
        
        sigma1=[c1(1,1)-c1(3,3) c2(1,1)-c2(3,3)];
        pa(k,2)=sum(sigma1);
        pa(k,1)=lamd;
        tt(i,2)=ac(i,2)+pa(k,2);
    else
        tt(i,2)=ac(i,2);
    end
    
    tt(i,1)=lamd;

  res(i)=(tt(i,2)-total(i,2))^2;

    
end
 resom=sum(res)/size(total,1)

plot(tt(:,1),tt(:,2))
plot(ac(:,1),ac(:,2))
plot(pa(:,1),pa(:,2))
