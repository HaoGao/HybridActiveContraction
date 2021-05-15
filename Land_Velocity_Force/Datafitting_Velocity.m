function f=Datafitting_Velocity(para)
global Velocity
%
r1=0.1475; r2=-0.3500; r3=1.1738;
lamd=2.3/1.85;
lmax=r1*lamd^2+r2*lamd+r3;%lmax=lmax*para(4);
ps=(1-lmax)*(9.25*lamd-1);

for i=1:size(Velocity,1)
    
     v=Velocity(i,1)*0.1*2.3/1.85;

     lamd=1.12;
     lmax1=r1*lamd^2+r2*lamd+r3;%lmax1=lmax1*para(4);
     lmax=lmax1*(1+para(2))/(1+para(2)*exp(-para(3)*v));

     po=((1-lmax)*(9.25*lamd-1))/ps;

     f(i)=po-Velocity(i,2);

end

end
