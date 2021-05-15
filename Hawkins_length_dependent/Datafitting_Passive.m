function f=Datafitting_Passive(para)
global passive
%%
%parapassive=[3.2790    0.0010    0.1418    4.6476];
%i=0; t0=0.01; tr=0.5;
%lmax=1.4;inte=0;lamd=lmax;pmax=104.72;k=0.25;vmax=193;kr=0.005;
 
 for i=1:size(passive,1)
     lamd=passive(i,1);
     [sigma_fs]=InteractionUM(lamd,para);

     f(i)=sigma_fs-passive(i,2);

 end
 
end
