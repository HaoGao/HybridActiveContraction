function dydt = oscfit(t,y)    

dydt = zeros(11,1);    % this creates an empty column     
%vector that you can fill with your two derivatives:   
%material parameters
a=y(4);b=y(5);c=y(6);
u1=y(7);u2=y(8);gama=y(9);
q=y(10);k=y(11);



dydt(1) = c*y(1)*(y(1)-a)*(1-y(1))-y(2)*y(1); 

phi=gama+u1*y(2)/(u2+y(1));

dydt(2) = phi*(-y(2)-c*y(1)*(y(1)-b-1));  

dydt(3) = q*y(1)-k*y(3);  

%In this case, y(1) is y1 and y(2) is y2, and dydt(1)     
%is dy1/dt and dydt(2) is dy2/dt. 
end
