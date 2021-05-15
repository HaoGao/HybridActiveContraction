function dydt = osc(t,y)    
dydt = zeros(3,1);    % this creates an empty column     
%vector that you can fill with your two derivatives:   
%material parameters
a=0.01;b=0.15;c=8;
u1=0.2;u2=0.3;gama=0.002;
q=0.1;k=q/1;

dydt(1) = c*y(1)*(y(1)-a)*(1-y(1))-y(2)*y(1); 

phi=gama+u1*y(2)/(u2+y(1));

dydt(2) = phi*(-y(2)-c*y(1)*(y(1)-b-1));  

dydt(3) = q*y(1)-k*y(3);  

%In this case, y(1) is y1 and y(2) is y2, and dydt(1)     
%is dy1/dt and dydt(2) is dy2/dt. 
end
