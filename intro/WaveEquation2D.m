%==========================================================================
% Finite element solution of wave equation using linear elements
% See also http://www.particleincell.com/blog/2012/matlab-fem/
% Meshgenerator, http://persson.berkeley.edu/distmesh
% D.E.: d^2u/dx^2 = 1/c^2 d^2u/dt^2
% B.C.: Dirichlet
%==========================================================================
clear all
close all
load meshdata1;
% c(x,y)===================================================================
c=@(X)1;
% =========================================================================
numberofnodes=size(p,1);
numberofelements=size(t,1);
% Initialisation of K, M ==================================================
K=zeros(numberofnodes,numberofnodes);
M=zeros(numberofnodes,numberofnodes);
% Computation of Ke, Me & assembly of K, M=================================
for e=1:numberofelements
   nodes=t(e,:); % row with nodenumbers of element e
   P=[ones(1,3);p(nodes,:)'];
   Coeffs=inv(P); 
   areaelemente=abs(det(P))/2; % area of element e
   diffPhi=Coeffs(:,2:end);
   Ke=diffPhi*diffPhi'*areaelemente;
   Me=1/3^2*ones(3)*1/(c(mean(p(nodes,:))))^2*areaelemente;
   K(nodes,nodes)=K(nodes,nodes)+Ke;
   M(nodes,nodes)=M(nodes,nodes)+Me;
end
Mlumped=diag(M(1:size(M,1)+1:end));
M=Mlumped;
% Dirichlet boundaries=====================================================
% =unique(boundedges(p,t));
Dirichlet=[
1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;23;24;33;40;49;56;65;72;81;88;97;103;
113;119;129;135;145;151;161;162;170;180;186;196;202;212;218;228;234;243;250;
259;266;275;282;291;298;307;308;315;316;317;318;319;320;321;322;323;324;325;
326;327;328;329;330];
% Put zeros and ones in Dirichlet DOFs=====================================
K(Dirichlet,:)=0;
K(Dirichlet,Dirichlet)=eye(numel(Dirichlet),numel(Dirichlet));
M(Dirichlet,:)=0;
M(Dirichlet,Dirichlet)=eye(numel(Dirichlet),numel(Dirichlet));
% Initial condition========================================================
C=zeros(numberofnodes,1);
centre1=find(((p(:,1)).^2+(p(:,2)-0.25).^2)<=0.20)';
x=p(centre1,1);
y=p(centre1,2)-0.25;
C(centre1)=exp(-20*(x.^2+y.^2));
Cdot=zeros(numberofnodes,1);
% Solving and plotting=====================================================
A=[zeros(numberofnodes),eye(numberofnodes);-M\K,zeros(numberofnodes)];
time=linspace(0,2,500);
xdot=@(T,Y)A*Y;
[time,Y]=ode45(xdot,time,[C;Cdot]);
Y=Y(:,1:numberofnodes); % only first half needed of [C;Cdot]===============
hold on;
view(150,20);
set(gcf,'color','w');
set(gca,'Projection','perspective');
camzoom(1.0);
axis([-1 2 -1 1 -0.5 1]);
axis equal;
axis vis3d;
axis on;
caxis([0 0.75]);
daspect([1 1 1/2])
for i=1:numel(time)
    cla;
    trisurf(t,p(:,1),p(:,2),Y(i,:),'facecolor','interp','Edgecolor',[0.5,0.5,0.5]);
    camlight right;
    lighting phong;
    drawnow; 
end
% Show node numbers========================================================
text(p(Dirichlet,1),p(Dirichlet,2),num2str(Dirichlet),'FontSize',9','Color',[0.75 0.75 0.75],...
'FontName','Cambria','HorizontalAlignment','left','VerticalAlignment','bottom');
rotate3d on;
