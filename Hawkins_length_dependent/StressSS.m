function [sigma1,sigma2]=StressSS(f0,n0,F,para)
global  record  paranum;
%dispersion
c1=[];c2=[];c3=[];c4=[];c5=[];c6=[];c7=[];c8=[];sigma1=[];sigma2=[];

%invariants
%f=[1 0 0]';s=[0 1 0]';
s0=[0 0 1]';

C=F'*F;
I_1=sum(diag(C));

I_4f=dot(f0,C*f0);
I_4_f=max(I_4f,1);

I_4n=dot(n0,C*n0);
I_4_n=max(I_4n,1);

I_4s=dot(s0,C*s0);
I_4_s=max(I_4s,1);

fxf=(F*f0)*(F*f0)';
nxn=(F*n0)*(F*n0)';
sxs=(F*s0)*(F*s0)';

%matrix
C=F'*F;
B=F*F';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%stress
a=para(1);b=para(2);
a_f=para(3);b_f=para(4);
%a_s=para(5);b_s=para(6);
%a_n=para(7);b_n=para(8);
%a_fs=para(9);b_fs=para(10);
%a_fn=para(11);b_fn=para(12);
%a_sn=para(13);b_sn=para(14);

%a_cfn=para(15);b_cfn=para(16);

%a=0;b=0;
%a_f=0;b_f=0;
a_s=0;b_s=0;
a_n=0;b_n=0;
a_fs=0;b_fs=0;
a_fn=0;b_fn=0;
a_sn=0;b_sn=0;
a_cfn=0;b_cfn=0;

    

    c1=(a*exp(b*(I_1-3)))*B;
    
    c2=(2*a_f*(I_4_f-1)*exp(b_f*(I_4_f-1)^2))*fxf;
  

    

    sigma1=[c1(1,1)-c1(3,3) c2(1,1)-c2(3,3)];
    sigma2=[c1(2,2)-c1(3,3) c2(2,2)-c2(3,3)];
    

 

 return
end
