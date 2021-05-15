
function [s11]=InteractionUM(lamd1,para)
global keyword record iter;

        thet=0;
        
        f0=[cos(thet) sin(thet) 0]';
        n0=[-sin(thet) cos(thet) 0]';
        F=zeros(3);
        k1=0;k2=0;
        
        lamd2=1/sqrt(lamd1);
        lamd3=1/sqrt(lamd1);
        F(1,1)=lamd1;F(2,2)=lamd2;F(3,3)=lamd3;

        [sigma1,sigma2]=StressSS(f0,n0,F,para);
           
        s11=sum(sigma1);
        %s22=s22/(n+1);

  

return

end

