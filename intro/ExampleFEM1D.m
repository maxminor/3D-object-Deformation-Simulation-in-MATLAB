%==========================================================================
% FEM solution of bvp: d^2/dx^2u + 4pi^2u = 0; u(0) = 0; du/dx(1) = 1 using
% quadratic 1D quadratic elements
% See also http://www.particleincell.com/blog/2012/finite-element-examples
%==========================================================================
clear all
close all
% Store coordinates in p===================================================
p=linspace(0,1,10)';
% Connectivity
t=[1 2
2 3
3 4
4 5
5 6
6 7
7 8
8 9
9 10];
%==========================================================================
TotalNumberOfNodes=size(p,1);
NumberOfElements=size(t,1);
% c(x), f(x) en lambda
c=@(x)1*ones(size(x));
f=@(x)0*ones(size(x));
labda=4*pi^2;
%========================================================================== 
% Quadratic elements contain 3 nodes per segment: add nodes in the middle
% of eacht segment
TotalNumberOfNodes=size(p,1);
NumberOfElements=size(t,1);
S=zeros(TotalNumberOfNodes);
counter=TotalNumberOfNodes+1;
for e=1:NumberOfElements
    nodes=t(e,:); 
    if (S(nodes(1),nodes(2))==0)
        S(nodes(1),nodes(2))=counter;
        S(nodes(2),nodes(1))=counter;
        p(counter,:)=mean(p([nodes(1) nodes(2)],:));
        counter=counter+1;
    end    
    t(e,3)=S(nodes(1),nodes(2));
end
% Update===================================================================
TotalNumberOfNodes=size(p,1);
NumberOfElements=size(t,1); 
% Initialisation of K, M and F=============================================
K=zeros(TotalNumberOfNodes); 
M=zeros(TotalNumberOfNodes); 
F=zeros(TotalNumberOfNodes,1); 
%==========================================================================
% lambda, lambda * two outer nodecoordinates yields integrationpoints
lambda=1/2*[
    1-1/sqrt(3) 1+1/sqrt(3);
    1+1/sqrt(3) 1-1/sqrt(3)];
% weights of Gaussian quadrature
w_Omega(1:2,1)=1;
hold on
% Loop over all elements
for e=1:NumberOfElements
    nodes=t(e,:); 
    % 3 by 3 matrix with rows: [ones;x;x^2]
    P=[ones(1,3);p(nodes,:)';p(nodes,1)'.^2];
    length=norm(diff(p(nodes([1:2]),:)));
    % Determine the two coordinates of the integration points on the 
    % edge of the element belonging to the Neuman boundary
    ip=lambda*p(nodes([1:2]),:);
    % plot integration points
    % plot(ip,0,'ko','MarkerSize',5,'MarkerFaceColor','y');
    % 3 by 2 matrix with rows: [ones;x;x^2] of two integration points
    I=[ones(1,2);ip';ip(:,1)'.^2];
    Phi=P\I;
    %========================================================================
    Ix=[zeros(1,2);ones(1,2);ip(:,1)'*2];
    diffI=Ix;
    diffPhi=P\diffI;
    %========================================================================
    cvalue=c(ip);
    fvalue=f(ip);
    Ke=zeros(size(P));
    Me=zeros(size(P));
    Fe_Omega=zeros(size(P,2),1);
    % Integrate and compute element stiffness matrix (3 by 3, 3 nodes)
    for i=1:size(Phi,2) 
        Ke=Ke+w_Omega(i)*cvalue(i)*diffPhi(:,i)*diffPhi(:,i)'*length/2;
        Me=Me+w_Omega(i)*Phi(:,i)*Phi(:,i)'*length/2;
        Fe_Omega=Fe_Omega+w_Omega(i)*fvalue(i)*Phi(:,i)*length/2;
    end
    % Assembly
    K(nodes,nodes)=K(nodes,nodes)+Ke;
    M(nodes,nodes)=M(nodes,nodes)+Me;
    F(nodes)=F(nodes)+Fe_Omega; 
end
% Proces Neumann boundary==================================================
F(10)=F(10)+c(1)*1;
% Proces Dirichlet boundary================================================
K(1,:)=0;
K(1,1)=1;
M(1,:)=0;
M(1,1)=1;
% Set the Dirichlet boundary (r.h.s.: value of fixed displament)===========
F(1)=0;
% Solve system=============================================================
U=(K-labda*M)\F;
% plot exact and FEM solution==============================================
box on
plot(p(1:10),1/2*sin(2*pi*p(1:10))/pi,'r','Linewidth',4)
plot(p(1:10),U(1:10),'o-g','Linewidth',2,'MarkerSize',8,'MarkerFaceColor','g');
