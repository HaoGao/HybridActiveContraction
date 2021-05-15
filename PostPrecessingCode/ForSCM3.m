%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
file='Job-Ca04s'
d0=LEstrain(file);

for cycle=1:12
    file=['cycle' num2str(cycle)]
     
	
    para=[0.25 0.1409 -0.3342 1.1210 0.0246 3.1338];

    %ComputeFFGFR;
    %the baseline values
  

    if cycle <=6
        dk=para(cycle)*0.1;
    else
        dk=para(cycle-6)*0.1;
    end
    d=[];
    d=LEstrain(file);
    dtest(:,cycle)=d;
    Smd(:,cycle)=(d-d0)/dk/2;
end
    
    
    
for i=1:6
    Sm=[];
    Sm=Smd(:,6+i)-Smd(:,i);
    Smbar(:,i)=Sm/norm(Sm);
    Smbv(:,i)=Sm;
    Smbvalue(i)=norm(Sm);
end

for i=1:6
    for j=1:6
        SCM(i,j)=Smbar(:,i)'*Smbar(:,j);
    end
end
        
