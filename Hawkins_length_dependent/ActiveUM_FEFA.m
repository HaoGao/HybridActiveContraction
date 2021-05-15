clear all; %close all;
global passive total tadd

load('EXPdataCau');
%load('lamda');

for ii=1:50

para0=1*rand(1,8);

%boundary values

Lb=[0.3, 0,  1e-3,0.0001,0.0001,0.0001,  0,0];
Ub=[1,   50, 50,50,50,50,   50,50];


options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_ActiveFEFA(para),para0,Lb,Ub,options);
para0=para;


[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_ActiveFEFA(para),para0,Lb,Ub,options);

abadd(ii)=resnorm;
dadd(ii,:)=para;

end

eges=find(abadd==min(abadd));
para=real(dadd(eges,:))

    a=para(3);b=para(4);
    a_f=para(5);b_f=para(6);
    
figure(1)
hold on

plot(total(:,1),total(:,2))
plot(active(:,1),active(:,2))
plot(passive(:,1),passive(:,2))

fc=1;
for i=1:size(total,1)
    lamd=total(i,1);
    
    F=[lamd 0 0; 0 1/sqrt(lamd) 0; 0 0 1/sqrt(lamd)];
    %lamd_a=call_osc4(lamd,para);
    %lamax=(para(1)*lamd+para(2));
    %lamax=(para(1)*lamd^2+para(2)*lamd+para(7));
    lamax=para(1);
     
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
    
    c1=(a*exp(b*(I_1-3)))*BE;
    c2=(2*a_f*(I_4_f-1)*exp(b_f*(I_4_f-1)^2))*fxfE;  

    sigma1=[c1(1,1)-c1(3,3) c2(1,1)-c2(3,3)];
    tt(i,2)=sum(sigma1)/detm;
    tt(i,1)=lamd;

  res(i)=(tt(i,2)-total(i,2))^2;

    
end
 resom=sum(res)/size(total,1)
 
 
k=0;
for i=1:size(total,1)
    lamd=total(i,1);

    F=[lamd 0 0; 0 1/sqrt(lamd) 0; 0 0 1/sqrt(lamd)];
%invariants
    f0=[1 0 0]';
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
        ac(i,2)=tt(i,2)-sum(sigma1);
    else
        ac(i,2)=tt(i,2);
    end
    ac(i,1)=lamd;
    
end

plot(tt(:,1),tt(:,2))
plot(ac(:,1),ac(:,2))
plot(pa(:,1),pa(:,2))
