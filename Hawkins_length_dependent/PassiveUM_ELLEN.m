clear all;%close all;
%
 global passive
 
 load('EXPdataCau');
 %dddd
for ii=1:1


para0=10*rand(1,5);

%boundary values
Lb=[0.01,0.1,0.01,0.01,0];
Ub=[10,10,20,10,1000];

options = optimoptions(@lsqnonlin,'Algorithm','trust-region-reflective');

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_Passive(para),para0,Lb,Ub,options);
para0=para;

[para,resnorm,residual]=lsqnonlin(@(para)Datafitting_Passive(para),para0,Lb,Ub,options);
para;


end
para

k=0;
for i=1:size(passive,1)
    lamd=passive(i,1);
        k=k+1;
        [sigma_fs]=InteractionUM(lamd,para);
        pa(k,1)=lamd;
        pa(k,2)=sigma_fs;
    
end



figure(3)
hold on

plot(passive(:,1),passive(:,2),'o')
plot(pa(:,1),pa(:,2))
    
