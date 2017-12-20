%==========================================================================
% FEM computation of deflection and stresses of a cantilever under uniform
% loads using quadratic elements
% see also http://www.particleincell.com/blog/2012/matlab-fem
%==========================================================================
clear all
close all
% Set tetrahedrical mesh===================================================
l=100;b=10;h=15;
X=linspace(0,b,3);
Y=linspace(0,l,20);
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
f=@(x,y,z)[0;0;-0.15];
% =========================================================================
TotalNumberOfNodes=size(p,1);
NumberOfElements=size(t,1);
% Linear elastic isotropic material========================================
D=E/((1+nu)*(1-2*nu))*[1-nu nu      nu      0 0 0
                       nu   1-nu    nu      0 0 0
                       nu   nu      1-nu    0 0 0
0 0 0 (1-2*nu)/2    0           0
0 0 0 0             (1-2*nu)/2  0
0 0 0 0             0           (1-2*nu)/2];
% Initialisation of K, M and F=============================================
K=zeros(3*TotalNumberOfNodes,3*TotalNumberOfNodes);
F=zeros(3*TotalNumberOfNodes,1);
%==========================================================================
for e=1:NumberOfElements
    NodesOfElement=t(e,:); % row containing nodenumbers of element e
    DegreesOfFreedom=reshape([3*NodesOfElement-2;3*NodesOfElement-1;3*NodesOfElement],1,3*numel(NodesOfElement));
    X=[ones(1,4);p(NodesOfElement,:)'];
    VolumeOfElement=abs(det(X))/6; % volume of element e
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
    % Computation of 12*12 Ke, 12*1 Fe
    Ke=B{e}'*D*B{e}*VolumeOfElement;
    Fe=1/4*repmat(eye(3),4,1)*f(mean(p(NodesOfElement,:)))*VolumeOfElement;
    K(DegreesOfFreedom,DegreesOfFreedom)=K(DegreesOfFreedom,DegreesOfFreedom)+Ke;
    F(DegreesOfFreedom)=F(DegreesOfFreedom)+Fe;
end
%==========================================================================
K(DegreesOfFreedomDirichlet,:)=0; % put zeros in boundary rows of K and F
K(DegreesOfFreedomDirichlet,DegreesOfFreedomDirichlet)=eye(numel(DegreesOfFreedomDirichlet));
% Set Dirichlet boundaries (value of fixed displacement)===================
F(DegreesOfFreedomDirichlet)=0;
% Solve system=============================================================
U=K\F;
Displacements=[U(1:3:end),U(2:3:end),U(3:3:end)];
Endpositions=p+magnification*Displacements;
% Computations of Nodestresses=============================================
MisesNodeStress=zeros(TotalNumberOfNodes,1);
for e=1:NumberOfElements
    NodesOfElement=t(e,:); % row containing nodenumbers of element e
    de=Displacements(t(e,:),:)';
    de=de(:); % [ux1,uy1,uz1,ux2,uy2,uz2,ux3,uy3,uz3,ux4,uy4,uz4]
    s=D*B{e}*de;
    mises=1./sqrt(2)*sqrt((s(1)-s(2)).^2+(s(2)-s(3)).^2+(s(3)-s(1)).^2+...
        6*s(4).^2+6*s(5).^2+6*s(6).^2);
    %mises=s(2);
    MisesNodeStress(NodesOfElement)=MisesNodeStress(NodesOfElement)+mises;
end
% Determine node frequencies (occurrences)=================================
NodeFrequencies=[];
for i=1:TotalNumberOfNodes
    NodeFrequencies=[NodeFrequencies;numel(find(t==i))];
end
MisesNodeStress=MisesNodeStress./NodeFrequencies;
%==========================================================================
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
        fill3(C2(:,1),C2(:,2),C2(:,3),MisesNodeStress(t(e,Faces(:,i))),'EdgeColor',[0,0,0]);
    end
end
% Label the nodes with a Dirichlet condition===============================
for i=1:numel(NodesDirichlet)
   quiver3(Endpositions(NodesDirichlet(i),1),Endpositions(NodesDirichlet(i),2),Endpositions(NodesDirichlet(i),3),5,0,0,'linewidth',1,'Color','m','MaxHeadSize',0.5);
   quiver3(Endpositions(NodesDirichlet(i),1),Endpositions(NodesDirichlet(i),2),Endpositions(NodesDirichlet(i),3),0,5,0,'linewidth',1,'Color','m','MaxHeadSize',0.5);
   quiver3(Endpositions(NodesDirichlet(i),1),Endpositions(NodesDirichlet(i),2),Endpositions(NodesDirichlet(i),3),0,0,5,'linewidth',1,'Color','m','MaxHeadSize',0.5);
end
axis equal
% Show node numbers========================================================
text(Endpositions(:,1),Endpositions(:,2),Endpositions(:,3),num2str([1:size(Endpositions,1)]'),'FontSize',8','Color',[0.5 0.5 0.5],...
'FontName','Cambria','HorizontalAlignment','left','VerticalAlignment','bottom');
rotate3d on;
%==========================================================================
MaximumDisplacement=min(Displacements(:,3));
% Compare to theoretical values
I=(b*h^3)/12;
w_q=f(x,y,z);
w_z=w_q(3)*b*h;
TheoreticalDeflectionUniformLoad=(w_z*l^4)/(8*E*I);
str=['FEM deflection: ',num2str(MaximumDisplacement),' vs. theoretical deflection:',num2str(TheoreticalDeflectionUniformLoad)];
disp(str)
