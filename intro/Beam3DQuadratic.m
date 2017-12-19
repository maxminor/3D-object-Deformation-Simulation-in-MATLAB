%==========================================================================
% FEM computation of deflection and stresses of a cantilever under uniform
% and point loads using quadratic elements
% See also http://www.particleincell.com/blog/2012/finite-element-examples
%==========================================================================
clear all
close all
% Set tetrahedrical mesh===================================================
l=100;b=10;h=15;
X=linspace(0,b,3);
Y=linspace(0,l,10);
Z=linspace(0,h,6);
[x,y,z]=meshgrid(X,Y,Z);
% Store coordinates in p===================================================
p=[x(:),y(:),z(:)];
t=delaunay(p); % connectivity
% Set parameters===========================================================
magnification=25;
E=2.0e5;
nu=0.3;
% List of nodes with a Dirichlet boundary condition (prescribed displacement)
NodesDirichlet=find(p(:,2)==0);
DegreesOfFreedomDirichlet=reshape([3*NodesDirichlet-2;3*NodesDirichlet-1;3*NodesDirichlet],1,3*numel(NodesDirichlet));
% f(x,y,z) force per unit volume ==========================================
f=@(x,y,z)[0;0;0];
% =========================================================================
TotalNumberOfNodes.old=size(p,1);
NumberOfElements=size(t,1);
S=zeros(TotalNumberOfNodes.old);
counter = TotalNumberOfNodes.old+1;
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
    if (S(nodes(1),nodes(4))==0)
        S(nodes(1),nodes(4))=counter;
        S(nodes(4),nodes(1))=counter;
        p(counter,:)=mean(p([nodes(1) nodes(4)],:));
        counter=counter+1;
    end
    if (S(nodes(2),nodes(4))==0)
        S(nodes(2),nodes(4))=counter;
        S(nodes(4),nodes(2))=counter;
        p(counter,:)=mean(p([nodes(2) nodes(4)],:));
        counter=counter+1;
    end
    if (S(nodes(3),nodes(4))==0)
        S(nodes(3),nodes(4))=counter;
        S(nodes(4),nodes(3))=counter;
        p(counter,:)=mean(p([nodes(3) nodes(4)],:));
        counter=counter+1;
    end
    t(e,5)=S(nodes(1),nodes(2));
    t(e,6)=S(nodes(2),nodes(3));
    t(e,7)=S(nodes(1),nodes(3));
    t(e,8)=S(nodes(1),nodes(4));
    t(e,9)=S(nodes(2),nodes(4));
    t(e,10)=S(nodes(3),nodes(4));
end
% Update===================================================================
TotalNumberOfNodes.new=size(p,1);
% Linear elastic isotropic material========================================
D=E/((1+nu)*(1-2*nu))*[1-nu nu      nu      0 0 0
                       nu   1-nu    nu      0 0 0
                       nu   nu      1-nu    0 0 0
0 0 0 (1-2*nu)/2    0           0
0 0 0 0             (1-2*nu)/2  0
0 0 0 0             0           (1-2*nu)/2];
% Initialisation of K and F================================================
K=zeros(3*TotalNumberOfNodes.new,3*TotalNumberOfNodes.new);
F=zeros(3*TotalNumberOfNodes.new,1);
%==========================================================================
% lambda, lambda * four corner nodecoordinates yields integrationpoints
lambda=[
    1/4+3/20*sqrt(5) 1/4-1/20*sqrt(5) 1/4-1/20*sqrt(5) 1/4-1/20*sqrt(5);
    1/4-1/20*sqrt(5) 1/4+3/20*sqrt(5) 1/4-1/20*sqrt(5) 1/4-1/20*sqrt(5);
    1/4-1/20*sqrt(5) 1/4-1/20*sqrt(5) 1/4+3/20*sqrt(5) 1/4-1/20*sqrt(5);
    1/4-1/20*sqrt(5) 1/4-1/20*sqrt(5) 1/4-1/20*sqrt(5) 1/4+3/20*sqrt(5)];
% weights of Gaussian quadrature
w_Omega(1:4,1)=1/4;
hold on
% Loop over all elements
for e=1:NumberOfElements
    nodes=t(e,:); % row containing nodenumbers of element e
    DegreesOfFreedom=reshape([3*nodes-2;3*nodes-1;3*nodes],1,3*numel(nodes));
    % 10 by 10 matrix with rows [1;x;y;z;|x^2;xy;xz;|y^2;yz;z^2] of the ten
    % nodes per element
    P=[ones(1,10);p(nodes,:)';...
        p(nodes,1)'.^2;p(nodes,1)'.*p(nodes,2)';p(nodes,1)'.*p(nodes,3)';...
        p(nodes,2)'.^2;p(nodes,2)'.*p(nodes,3)';p(nodes,3)'.^2];
    volume=abs(det(P(1:4,1:4)))/6;
    % Integration points
    ip=lambda*p(nodes(1:4),:); 
    % plot the integration points
    % plot3(ip(:,1),ip(:,2),ip(:,3),'ro','MarkerSize',2,'MarkerFaceColor','r');
    
    % I is a 10 by 1 matrix with rows [ones;x;y;z;|x^2;xy;xz;|y^2;yz;z^2] of
    % the four integration points per element
    I=[ones(1,4);ip';...
        ip(:,1)'.^2;ip(:,1)'.*ip(:,2)';ip(:,1)'.*ip(:,3)';...
        ip(:,2)'.^2;ip(:,2)'.*ip(:,3)';ip(:,3)'.^2]; 
    Phi=P\I; % Phi(ipx,ipy,ipz)
    Ix=[zeros(1,4);ones(1,4);zeros(1,4);zeros(1,4);...
        2*ip(:,1)';ip(:,2)';ip(:,3)';
        zeros(1,4);zeros(1,4);zeros(1,4)];
    Iy=[zeros(1,4);zeros(1,4);ones(1,4);zeros(1,4);...
        zeros(1,4);ip(:,1)';zeros(1,4);
        2*ip(:,2)';ip(:,3)';zeros(1,4)];
    Iz=[zeros(1,4);zeros(1,4);zeros(1,4);ones(1,4);...
        zeros(1,4);zeros(1,4);ip(:,1)';
        zeros(1,4);ip(:,2)';2*ip(:,3)'];
    diffI=[Ix Iy Iz];
    diffPhi=P\diffI; % diffPhi(ipx,ipy,ipz)
    %========================================================================
    Ke=zeros(3*size(P));
    Fe_Omega=zeros(3*size(P,2),1);
    fvalue=f(ip(:,1),ip(:,2),ip(:,3));
    for j=1:4 % loop over integration points       
        B{e}=[];
        for i=1:10 % loop over 10 element nodes
            Bsmall=[
                diffPhi(i,j)    0                               0;
                0                       diffPhi(i,j+4)          0;
                0                       0                       diffPhi(i,j+8);
                diffPhi(i,j+4)          diffPhi(i,j)            0;
                0                       diffPhi(i,j+8)          diffPhi(i,j+4);
                diffPhi(i,j+8)          0                       diffPhi(i,j)];
            B{e}=[B{e},Bsmall];
        end
        % Integrate and compute element stiffness matrix (30 by 30, 10 nodes times three DOF)
        Ke=Ke+w_Omega(j)*B{e}'*D*B{e}*volume;
        %==================================================================
        FE_OMEGA=[];
        for k=1:numel(Phi(:,j))
            % Integrate and compute element load vector (30 by 1, 10 nodes times three DOF)
            FE_OMEGA=[FE_OMEGA;w_Omega(j)*fvalue*Phi(k,j)*volume];
        end
        Fe_Omega=Fe_Omega+FE_OMEGA;
    end
    % Assembly
    K(DegreesOfFreedom,DegreesOfFreedom)=K(DegreesOfFreedom,DegreesOfFreedom)+Ke;
    F(DegreesOfFreedom)=F(DegreesOfFreedom)+Fe_Omega;
end
%==========================================================================
K(DegreesOfFreedomDirichlet,:)=0; % put zeros in boundary rows of K and F
K(DegreesOfFreedomDirichlet,DegreesOfFreedomDirichlet)=eye(numel(DegreesOfFreedomDirichlet));
% Set Dirichlet boundaries (value of fixed displacement)===================
F(DegreesOfFreedomDirichlet)=0;
%==========================================================================
% Apply a pointload at specified nodes
NodesWithPointLoad=[30,10,180,160];
PointLoads=[[0;0;-500/4],[0;0;-500/4],[0;0;-500/4],[0;0;-500/4]];
F([3*NodesWithPointLoad-2;3*NodesWithPointLoad-1;3*NodesWithPointLoad])=PointLoads;
% Solve system=============================================================
U=K\F;
%==========================================================================
Displacements=[U(1:3:end),U(2:3:end),U(3:3:end)];
% Compute displacement of corners of element
Displacements=Displacements(1:TotalNumberOfNodes.old,:);
Endpositions=[x(:),y(:),z(:)]+magnification*Displacements;
% Computations of Nodestresses (Constant strain!)==========================
MisesNodeStress=zeros(TotalNumberOfNodes.old,1);
for e=1:NumberOfElements
    nodesofcorners=t(e,1:4); % row containing four nodenumbers of corners of element e
    X=[ones(1,4);p(nodesofcorners,:)'];
    Coeffs=inv(X);
    diffPhi=Coeffs(:,2:end);
    B{e}=[];
    for i=1:4
        Bsmall=[
                diffPhi(i,1)    0               0;
                0               diffPhi(i,2)    0;
                0               0               diffPhi(i,3);
                diffPhi(i,2)    diffPhi(i,1)    0;
                0               diffPhi(i,3)    diffPhi(i,2);
                diffPhi(i,3)    0               diffPhi(i,1)];
        B{e}=[B{e},Bsmall];
    end
    de=Displacements(nodesofcorners,:)';
    de=de(:); % [ux1;uy1;uz1;ux2;uy2;uz2;ux3;uy3;uz3;ux4;uy4;uz4]
    s=D*B{e}*de;
    mises=1./sqrt(2)*sqrt((s(1)-s(2)).^2+(s(2)-s(3)).^2+(s(3)-s(1)).^2+...
        6*s(4).^2+6*s(5).^2+6*s(6).^2);
    %mises=s(2);
    MisesNodeStress(nodesofcorners)=MisesNodeStress(nodesofcorners)+mises;
end
% Determine corner node frequencies (occurrences)==========================
NodeFrequencies=[];
for i=1:TotalNumberOfNodes.old
    NodeFrequencies=[NodeFrequencies;numel(find(t(:,1:4)==i))];
end
MisesNodeStress=MisesNodeStress./NodeFrequencies;
% Plotting=================================================================
set(gcf,'color','w')
set(gca,'Projection','perspective')
view(-220,5)
axis vis3d
axis on
hold on
xlabel('x')
ylabel('y')
zlabel('z')
handle=colorbar;
ylabel(handle,'Von Mises \sigma')
%==========================================================================
Faces=[1 1 2 3;2 2 3 1;3 4 4 4;1 1 2 3];
for e=1:NumberOfElements
    for i=1:size(Faces,1)
        C1=p(t(e,Faces(:,i)),:);
        C2=Endpositions(t(e,Faces(:,i)),:);
        % Show original mesh
        fill3(C1(:,1),C1(:,2),C1(:,3),0,'LineStyle','--','EdgeColor',[0.5,0.5,0.5],'FaceColor','none');
        % Show deformed mesh
        fill3(C2(:,1),C2(:,2),C2(:,3),MisesNodeStress(t(e,Faces(:,i))),'EdgeColor',[0,0,0],'facealpha',1);
    end
end
% Label the nodes with a Dirichlet condition===============================
for i=1:numel(NodesDirichlet)
   quiver3(Endpositions(NodesDirichlet(i),1),Endpositions(NodesDirichlet(i),2),Endpositions(NodesDirichlet(i),3),5,0,0,'linewidth',1,'Color','m','MaxHeadSize',0.5);
   quiver3(Endpositions(NodesDirichlet(i),1),Endpositions(NodesDirichlet(i),2),Endpositions(NodesDirichlet(i),3),0,5,0,'linewidth',1,'Color','m','MaxHeadSize',0.5);
   quiver3(Endpositions(NodesDirichlet(i),1),Endpositions(NodesDirichlet(i),2),Endpositions(NodesDirichlet(i),3),0,0,5,'linewidth',1,'Color','m','MaxHeadSize',0.5);
end
% Plot pointloads
scalefactor=0.05;
for i=1:numel(NodesWithPointLoad)
   quiver3(Endpositions(NodesWithPointLoad(i),1),Endpositions(NodesWithPointLoad(i),2),Endpositions(NodesWithPointLoad(i),3),...
   scalefactor*PointLoads(1,i),scalefactor*PointLoads(2,i),scalefactor*PointLoads(3,i),...
   'linewidth',2,'Color',[1 0.5 0],'MaxHeadSize',3);
end
axis equal
% Show node numbers========================================================
text(Endpositions(:,1),Endpositions(:,2),Endpositions(:,3),num2str([1:size(Endpositions,1)]'),'FontSize',10','Color',[0.75 0.75 0.75],...
'FontName','Cambria','HorizontalAlignment','left','VerticalAlignment','bottom');
rotate3d on;
%==========================================================================0
MaximumDisplacement=max(abs(Displacements(:,3)));
% Compare to theoretical values
I=(b*h^3)/12;
w_q=f(x,y,z);
w_z=w_q(3)*b*h;
TheoreticalDeflectionUniformLoad=(w_z*l^4)/(8*E*I);
w_P=-500;
TheoreticalDeflectionPointLoad=abs(w_P*l^3)/(3*E*I);
str=['FEM deflection: ',num2str(MaximumDisplacement),' vs. theoretical deflection:',num2str(TheoreticalDeflectionUniformLoad+TheoreticalDeflectionPointLoad)];
disp(str)

