clear all; 
%close all;
global Velocity

load('Land_Velocity2');
expds=load('Land_VelocityEXP');
%Velocity=expds.Velocity;

 
for st=1:10
para0=10*rand(1,9);

%boundary values 0.7559
Lb=[0, 0,   0,   0.05,  1, 1,  0, 0, 0];
Ub=[1, 500,  200,  1.5, 2, 50, 1, 10, 2];

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_Velocity(para),para0,Lb,Ub,options);
para0=para;

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_Velocity(para),para0,Lb,Ub,options);
para;

sdsf(st)=resnorm;
adasd(st,:)=para;
end
min(sdsf)
find(sdsf==ans)
para=adasd(ans,:)



lamd=2.3/1.85;
r1=0.1475; r2=-0.3500; r3=1.1738;

fc=1;
lmax=r1*lamd^2+r2*lamd+r3;%lmax=lmax*para(4);
ps=(1-lmax)*(9.25*lamd-1);


for i=1:size(Velocity,1)

     v=Velocity(i,1)*0.1*2.3/1.85;

     lamd=1.12;
     lmax1=r1*lamd^2+r2*lamd+r3; 
     %lmax1=lmax1*para(4);
     lmax=lmax1*(1+para(2))/(1+para(2)*exp(-para(3)*v));

   
     po=((1-lmax)*(9.25*lamd-1))/ps;
   
     ve(i,2)=po;
     ve(i,1)=Velocity(i,1)*0.23/1.85;
     %end
     
 end

 figure(6)
 hold on
 plot(ve(:,1),ve(:,2))
 plot(Velocity(:,1)*0.23/1.85,Velocity(:,2),'O')
 plot(expds.Velocity(:,1)*0.23/1.85,expds.Velocity(:,2),'*')
 
 para(2)=0.02; para(3)=3;
 i=0;
 for ii=0:1000
     i=i+1;
     v=ii*0.01;
     lamd=1.12;
     lmax1=r1*lamd^2+r2*lamd+r3; %lmax1=lmax1*para(4);
     
     %lmax=lmax1*(1+para(2)*log(v*para(3)+1));
     %lmax=lmax1*(1+para(2)/para(3)-para(2)/(v+para(3)));
    %lmax=lmax1*(1+para(2)*(1-exp(-para(3)*v))/(1+exp(-para(3)*v))/2);
    
     lmax=lmax1*(1+para(2))/(1+para(2)*exp(-para(3)*v));
     po=((1-lmax)*(9.25*lamd-1))/ps;

     %lamd_a(i)=real(xi*lmax/(1+fc*(xi-1)));
     lamd_a(i)=(1+para(2))/(1+para(2)*exp(-para(3)*v));
   
     ved(i)=po;
     gdv(i)=v;
 end
figure(7)
hold on
 plot(lamd_a)
 plot(ved)
 figure
 plot(gdv,lamd_a)
 
% load('GU-Force');
% plot(force(:,1)-50,force(:,2)*50)
     