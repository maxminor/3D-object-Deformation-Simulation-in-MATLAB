%==========================================================================
% Finite element solution of 2D heat equation (+ animation) using quadratic elements
% See also http://www.particleincell.com/blog/2012/matlab-fem
% Meshgenerator, http://persson.berkeley.edu/distmesh
% D.E.: d^2T/dx^2 + d^2T/dy^2 + f(x,y) = dT/dt
% B.C.: Neumann: dT/dn (Neumann) = -1
% B.C.: Dirichlet: corner of omega: 1, rest 0
%==========================================================================
clear all
close all
clc
%==========================================================================
load meshdata2
p(1,:) =[0,0.25];
p(18,:)=[0.25,0];
p(54,:)=[1,0.25];
p(57,:)=[1,0.75];
% Neumann boundary=========================================================
Neumann_e=[
    57 55;
    55 58;
    58 54];
% Dirichlet boundary=======================================================
Dirichlet_e=[
    1 8;
    8 5;
    5 12;
    12 18;
    4 2;
    2 6;
    6 3;
    3 7;
    7 11;
    11 17;
    17 25;
    25 33;
    33 40;
    40 46;
    46 53;
    53 59;
    59 57;
    54 56;
    56 52;
    52 49;
    49 42;
    42 34;
    34 26];
%========================================================================== 
c=@(x,y) 1*ones(size(x));
f=@(x,y) 10*sin(pi*x).*sin(pi*y);
%========================================================================== 
% Quadratic elements contain six nodes per triangle: add three nodes at 
% the middle of the edges of each triangle
TotalNumberOfNodes.old=size(p,1);
NumberOfElements=size(t,1);
S=zeros(TotalNumberOfNodes.old);
counter=TotalNumberOfNodes.old+1;
for e=1:NumberOfElements
    nodes=t(e,:); 
    if (S(nodes(1),nodes(2))==0)
        S(nodes(1),nodes(2))=counter;
        S(nodes(2),nodes(1))=counter;
        p(counter,:)=mean(p([nodes(1) nodes(2)],:));
        counter=counter+1;
    end
    if (S(nodes(2),nodes(3))==0)
        S(nodes(2),nodes(3))=counter;
        S(nodes(3),nodes(2))=counter;
        p(counter,:)=mean(p([nodes(2) nodes(3)],:));
        counter=counter+1;
    end
    if (S(nodes(1),nodes(3))==0)
        S(nodes(1),nodes(3))=counter;
        S(nodes(3),nodes(1))=counter;
        p(counter,:)=mean(p([nodes(1) nodes(3)],:));
        counter=counter+1;
    end
    t(e,4)=S(nodes(1),nodes(2));
    t(e,5)=S(nodes(2),nodes(3));
    t(e,6)=S(nodes(1),nodes(3));
end
% Add midpoints to the edges of the Dirichlet and Neumann boundaries=======
Dirichlet=[];
for i=1:size(Dirichlet_e)
    nodes=Dirichlet_e(i,:);
    Dirichlet=[Dirichlet;[nodes(1),S(nodes(1),nodes(2)),nodes(2)]];
end
Neumann=[];
for i=1:size(Neumann_e)
    nodes=Neumann_e(i,:);
    Neumann=[Neumann;[nodes(1),S(nodes(1),nodes(2)),nodes(2)]];
end
% Update===================================================================
TotalNumberOfNodes.new=size(p,1);
% Initialisation of K, M and F==============================================
K=zeros(TotalNumberOfNodes.new);
M=zeros(TotalNumberOfNodes.new); 
F=zeros(TotalNumberOfNodes.new,1); 
%========================================================================== 
% lambda, lambda * three corner nodecoordinates yields integration points
lambda=[
    2/3 1/6 1/6;
    1/6 2/3 1/6;
    1/6 1/6 2/3];
w_Omega(1:3,1)=1/3;
hold on
for e=1:NumberOfElements
    nodes=t(e,:); 
    % 6 by 6 matrix with rows: [ones;x;y;x^2;xy;y^2]:
    P=[ones(1,6);p(nodes,:)';p(nodes,1)'.^2;p(nodes,1)'.*p(nodes,2)';p(nodes,2)'.^2];
    areaelemente=abs(det(P(1:3,1:3)))/2; 
    % three integration points within each triangle
    ip=lambda*p(nodes(1:3),:); 
    % plot the integration points
    plot(ip(:,1),ip(:,2),'ro','MarkerSize',4,'MarkerFaceColor','r');
    %========================================================================
    % 6 by 3 matrix with rows: [ones;x;y;x^2;xy;y^2]' (of three integration
    % points)
    I=[ones(1,3);ip';ip(:,1)'.^2;ip(:,1)'.*ip(:,2)';ip(:,2)'.^2]; 
    Phi=P\I; % Phi(ip1,ip2,ip3)
    %========================================================================
    Ix=[zeros(1,3);ones(1,3);zeros(1,3);ip(:,1)'*2;ip(:,2)';zeros(1,3)];
    Iy=[zeros(1,3);zeros(1,3);ones(1,3);zeros(1,3);ip(:,1)';ip(:,2)'*2];
    diffI=[Ix Iy];
    diffPhi=P\diffI;
    %========================================================================
    cvalue=c(ip(:,1),ip(:,2));
    fvalue=f(ip(:,1),ip(:,2));
    Ke=zeros(size(P));
    Me=zeros(size(P));
    Fe_Omega=zeros(size(P,2),1);
    % Computation of element matrices
    for i=1:size(diffPhi,2)
        % modulo because of d/dx, and d/dy, see also diffI (w_Omega needed two laps)
        Ke=Ke+w_Omega(mod(i-1,3)+1)*cvalue(mod(i-1,3)+1)*diffPhi(:,i)*diffPhi(:,i)'*areaelemente;
    end
    for i=1:size(Phi,2)
        Me=Me+w_Omega(i)*Phi(:,i)*Phi(:,i)'*areaelemente;
        %Fe_Omega=Fe_Omega+w_Omega(i)*fvalue(i)*Phi(:,i)*areaelemente;
    end
    Fe_Omega= Phi*(w_Omega.*fvalue)*areaelemente;
    % Assembly
    K(nodes,nodes)=K(nodes,nodes)+Ke;
    M(nodes,nodes)=M(nodes,nodes)+Me;
    F(nodes)=F(nodes)+Fe_Omega;
end  
K(Dirichlet,:)=0; 
K(Dirichlet,Dirichlet)=eye(numel(Dirichlet));
M(Dirichlet,:)=0;
M(Dirichlet,Dirichlet)=eye(numel(Dirichlet));
% Proces Neumann boundary==================================================
% lambda, lambda * two corner nodecoordinates yields integrationpoints
lambda=1/2*[
    1-1/sqrt(3) 1+1/sqrt(3);
    1+1/sqrt(3) 1-1/sqrt(3)];
% weights of Gaussian quadrature
w_Gamma(1:2,1)=1;
for e=1:size(Neumann,1)
   nodes=Neumann(e,:);
   % 6 by 3 matrix with rows: [1 x y x^2 xy y^2]' (three nodes on 
   % the edge of each element)
   P=[ones(1,3);p(nodes,:)';p(nodes,1)'.^2;p(nodes,1)'.*p(nodes,2)';p(nodes,2)'.^2];
   lengte=norm(diff(p(nodes([1,3]),:)));
   % Determine the two coordinates of the integration points on the edge
   % of the element belonging tot the Neumann boundary
   ip=lambda*p(nodes([1,3]),:); 
   plot(ip(:,1),ip(:,2),'ko','MarkerSize',5,'MarkerFaceColor','y');
   % 6 by 2 matrix with rows: [ones;x;y;x^2;xy;y^2] by the
   % two integration points on the Neumann element edge
   I=[ones(1,2);ip';ip(:,1)'.^2;ip(:,1)'.*ip(:,2)';ip(:,2)'.^2];
   Phi=P\I;
   % Computation of element matrices
   Fe_Gamma=zeros(size(P,2),1);
   cvalue=c(ip(:,1),ip(:,2));
   for i=1:size(Phi,2)        
        %(-1), gradient of temperature field at Neumann boundary
        Fe_Gamma=Fe_Gamma+w_Gamma(i)*cvalue(i)* (-1) *Phi(:,i)*lengte/2;
   end
   F(nodes)=F(nodes)+Fe_Gamma;
end
F(Dirichlet)=0;
F(Dirichlet(1:4,:))=1;
T=K\F; 
% Plot the solution========================================================
set(gcf,'color','w')
set(gca,'Projection','perspective')
trisurf(t(:,1:3),p(:,1),p(:,2),0*T,'facecolor','none');
trisurf(t(:,1:3),p(:,1),p(:,2),T,'facecolor','interp');
text(p(:,1),p(:,2),T,num2str([1:size(p,1)]'),'FontSize',10','Color',[0.5,0.5,0.5],...
'FontName','Cambria','HorizontalAlignment','left','VerticalAlignment','bottom');
caxis([0 1])
daspect([1 1 2])
view(70,20)
axis vis3d
axis([0 1 0 1 0 1]);
axis on;
drawnow;
%==========================================================================
% Initial condition========================================================
Tstart=zeros(TotalNumberOfNodes.new,1);
Tstart(Dirichlet(1:4,:))=1;
%==========================================================================
% Plot the solution========================================================
figure(2)
set(gcf,'color','w');
set(gca,'Projection','perspective');
caxis([0 1]);
daspect([1 1 2]);
view(70,20);
axis vis3d;
axis([0 1 0 1 0 1]);
axis on;
hold on;
time=linspace(0,0.25,20);
xdot=@(time,T) M\(F-K*T);
[~,T]=ode23(xdot,time,Tstart);
hold on;
zlabel('T(x,y,t)');
for i=1:numel(time)
    cla;
    trisurf(t(:,1:3),p(:,1),p(:,2),T(i,:),'facecolor','interp');
    drawnow;
    %mov(i) = getframe(gcf);
end
%movie2avi(mov, '2D Heating.avi', 'compression', 'None');