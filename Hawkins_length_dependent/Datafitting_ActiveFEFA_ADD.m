function f=Datafitting_ActiveFEFA_ADD(para)
global passive total active
%%
%parapassive=[3.2790    0.0010    0.1418    4.6476];
%i=0; t0=0.01; tr=0.5;
%lmax=1.4;inte=0;lamd=lmax;pmax=104.72;k=0.25;vmax=193;kr=0.005;
k=0;
for i=1:size(passive,1)
    lamd=passive(i,1);
    
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
    a=para(3);b=para(4);
    a_f=para(5);b_f=para(6);
    
    if I_4_f>1
        
        k=k+1;
        c1=(a*exp(b*(I_1-3)))*B;
        c2=(2*a_f*(I_4_f-1)*exp(b_f*(I_4_f-1)^2))*fxf;
        sigma1=[c1(1,1)-c1(3,3) c2(1,1)-c2(3,3)];
        f(k)=sum(sigma1)-passive(i,2);
    end

    
end

 fc=1;
for i=1:size(active,1)
    k=k+1;
    lamd=active(i,1);
    
    F=[lamd 0 0; 0 1/sqrt(lamd) 0; 0 0 1/sqrt(lamd)];
    %ta=call_osc4(lamd,para);
    %lamd_a=call_osc4(lamd,para);
    %lamax=(para(1)*lamd+para(2));
    %lamax=(para(8)*lamd^2+para(9)*lamd+para(10));
    lamax=para(8);
     
     ps=(0-1)/(0-lamax);
     %ta=lamax/((1-fc)*lamax+fc);
    lamd_a=ps*lamax/(1+fc*(ps-1));

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
    
    eps=para(7);
   
    c2=2*eps*(I_4_f-1)*fxfE;  
    sigma1=[];
    sigma1=[c2(1,1)-c2(3,3)];
    f(k)=sum(sigma1)/detm-active(i,2);
    
end

end
