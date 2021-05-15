clear all; %close all;

file='Job-Ca04s_NoV';
lamdinter=load('D:\DBGuan\ActiveStrian_Github\ActiveStrain\ABAQUS\fibreshortening3.txt');
node=load(['D:\DBGuan\ActiveStrian_Github\ActiveStrain\ABAQUS\' file 'DFALL.txt']);
lamd=load(['D:\DBGuan\ActiveStrian_Github\ActiveStrain\ABAQUS\' file 'MeanLamd.txt']);
apex=load('D:\DBGuan\ActiveStrian_Github\ActiveStrain\ABAQUS\topregion.txt');
element_d=load('Edata.txt');
node0=load('Ndata.txt');
fs_di=load('fibresheet.txt');
fs_dir=fs_di(:,2:7); 
%%***********************************************************************************************************
%Compute total F
%read the message in the input file

node_d=node(1:26010,:);
% re order
idx=[];
[~,idx] = sort(node_d(:,1)); % sort just the first column
node_d = node_d(idx,:);
node_dxdydz_d=node_d-node0;

for i=1:size(element_d,1)
        xyztet=[];
        for j=1:4
            xyztet(1,j)=node0(element_d(i,j+1),1+1);
            xyztet(2,j)=node0(element_d(i,j+1),2+1);
            xyztet(3,j)=node0(element_d(i,j+1),3+1);
        
            dxdydz(1,j)=node_dxdydz_d(element_d(i,j+1),1+1);
            dxdydz(2,j)=node_dxdydz_d(element_d(i,j+1),2+1);
            dxdydz(3,j)=node_dxdydz_d(element_d(i,j+1),3+1);
        end
    
        [abc, Vcol]=IsoTet4ShapeFunDer(xyztet);
    
	%the deformation gradient tensor
        F=[];
        F=eye(3)+(dxdydz*abc)/6/Vcol; % defined in global cardisen coordinate
	% determinant of F
		F=F*(det(F))^(-1/3);
        s0=fs_dir(i,4:6)';
        f0=fs_dir(i,1:3)';
        FS(i,:)=F*s0/norm(F*s0);
        FF(i,:)=F*f0/norm(F*f0);
        
end


for frame=1:21
    
    num2str(frame)
    
    k1=(frame-1)*26010+1; k2=frame*26010;
    node_dxdydz_d=[];
    node_dxdydz_d=node(k1:k2,:)-node(1:26010,:);
    node_dxdydz_d(:,1)=node(1:26010,1);
    idx=[];
    [~,idx] = sort(node_dxdydz_d(:,1)); % sort just the first column
    node_dxdydz_d = node_dxdydz_d(idx,:);
    ecc=0;err=0;ell=0;
    for i=1:size(element_d,1)
        if ismember(i,apex)
        
        xyztet=[];
        for j=1:4
            xyztet(1,j)=node_d(element_d(i,j+1),1+1);
            xyztet(2,j)=node_d(element_d(i,j+1),2+1);
            xyztet(3,j)=node_d(element_d(i,j+1),3+1);
        
            dxdydz(1,j)=node_dxdydz_d(element_d(i,j+1),1+1);
            dxdydz(2,j)=node_dxdydz_d(element_d(i,j+1),2+1);
            dxdydz(3,j)=node_dxdydz_d(element_d(i,j+1),3+1);
        end
    
        [abc, Vcol]=IsoTet4ShapeFunDer(xyztet);
    
	%the deformation gradient tensor
        F=[];
        F=eye(3)+(dxdydz*abc)/6/Vcol; % defined in global cardisen coordinate
	% determinant of F
		F=F*(det(F))^(-1/3);

        lon=[0 0 1]';
        %r=F*FS(i,:)';
        r=FS(i,:)';
        r=r/norm(r);
        %cf1=F*FF(i,:)';
        cf1=FF(i,:)';
        cf=[cf1(1) cf1(2) 0]';
        c=cf/norm(cf);

        %c=cross(lon,r);
        %c=c/norm(c);
        %lon=cross(r,c);
        %lon=lon/norm(lon);
        
        r=cross(c,lon);
        r=r/norm(r);
        
        C=F'*F;
        E=(C-eye(3))/2;
        ec=dot(c,E*c);
        er=dot(r,E*r);
        el=dot(lon,E*lon);
        
        ecc=ecc+ec;
        err=err+er;
        ell=ell+el;
       
        end
    end
    nor=size(apex,1)*size(apex,2);%size(element_d,1)-
    ecrl(frame,1)=ecc/nor;
    ecrl(frame,2)=err/nor;
    ecrl(frame,3)=ell/nor;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%    fid2 = fopen([file 'DF_CRL.txt'],'w');
%    for i = 1 : size(ecrl,1)
%        fprintf(fid2, '\t%14.10f,\t%14.10f,\t%14.10f\n', ...
%            ecrl(i,1), ecrl(i,2), ecrl(i,3));
%    end
%    fclose(fid2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:size(lamd,1)
    b1(i)=(lamd(i,2)-lamd(1,2))/lamd(1,2)*100;
end

figure(1)
hold on
%plot(ecrl(:,1),b1,'o')
%plot(ecrl(:,2),b1,'o')
%plot(ecrl(:,3),b1,'o')

%%%%%%%%%%%%%%%%%%%%%%%

E = 0:0.1:2;
xi=0:0.01:2;
yispl=interp1(E,ecrl(:,1),xi,'spline');
yisp2=interp1(E,ecrl(:,2),xi,'spline');
yisp3=interp1(E,ecrl(:,3),xi,'spline');
yispb=interp1(E,b1,xi,'spline');


plot(yispl(1:8),yispb(1:8), '--',yispl(8:30),yispb(8:30), ':',yispl(30:201),yispb(30:201), '-')
plot(yisp2(1:8),yispb(1:8), '--',yisp2(8:30),yispb(8:30), ':',yisp2(30:201),yispb(30:201), '-')
plot(yisp3(1:8),yispb(1:8), '--',yisp3(8:30),yispb(8:30), ':',yisp3(30:201),yispb(30:201), '-')

    
plot(ecrl(6,1),b1(6),'o')
plot(ecrl(6,2),b1(6),'o')
plot(ecrl(6,3),b1(6),'o')