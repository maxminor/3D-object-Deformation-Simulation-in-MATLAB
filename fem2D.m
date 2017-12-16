function fem2D()

%create square mesh
[tri, v] = createSquareMesh(10,10, 1,1);
disp(tri);
% [tri, v] = createCircleMesh(5,0.01);
trimesh(tri, v(:,1), v(:,2));
numVerts = size(v,1);
pos = zeros(numVerts,2);
ylim([-5 5]);
xlim([-10 10]);
%floor is at zero
odeFun = @(time, state)(femOde(time,state,tri,v,pos));
outputFun = @(time, state, flag)(femOutputFcn(time,state,flag,pos, tri, v));
options = odeset('OutputFcn', outputFun);

y0 = [reshape(v', 2*numVerts,1); zeros(2*numVerts,1)];
[time, state] = ode45(odeFun, [0, 10], y0, options);

pos(1:numVerts,:) = [state(end, 1:2:2*numVerts)' state(end, 2:2:2*numVerts)'];
trimesh(tri, pos(:,1), pos(:,2));
    
end

function y = femOde(time, state, tri, v, pos)
    numVerts = size(v,1);
    pos(1:numVerts,:) = [state(1:2:2*numVerts) state(2:2:2*numVerts)];
    
    %hack to fix object to the floor
    %just set y positions of vertices to zero
    floorVert = abs(v(:,2)) < 1e-6; %some arbitrary toleranc
    
    %set positions to zero
    pos(floorVert,2) = 0;
    
    %set velocities to zero
    vel = [state((2*numVerts+1):2:end) state((2*numVerts+2):2:end)];
    vel(floorVert,2) = 0;
    y = [reshape(vel', 2*numVerts,1) ; femAccelerations(tri, v, pos)];
    
end

%compute FEM accelerations
function a = femAccelerations(tri, v, pos)

numTris = size(tri, 1);
numVerts = size(v,1);
f = zeros(numVerts,2);
m = zeros(numVerts, 2);
density = 1.0;
g = [0 -9.8]';

for i=1:numTris
    
    %Edge vectors
    E1 = (v(tri(i,2), :)-v(tri(i,1), :))';
    E2 = (v(tri(i,3), :)-v(tri(i,1), :))';
    
    e1 = (pos(tri(i,2), :)-pos(tri(i,1), :))';
    e2 = (pos(tri(i,3), :)-pos(tri(i,1), :))';
    e3 = (pos(tri(i,3), :)-pos(tri(i,2), :))';
    area = 0.5*det([e1 e2]);
    
    me = density*area;
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %compute F
%     disp(pos);
%     cref = [(v(tri(i,1), 1)+v(tri(i,2), 1) + v(tri(i,3), 1))/3 , (v(tri(i,1), 2)+v(tri(i,2), 2) + v(tri(i,3), 2))/3 ]';
%     cpos = [(pos(tri(i,1), 1)+pos(tri(i,2), 1) + pos(tri(i,3), 1))/3 , (pos(tri(i,1), 2)+pos(tri(i,2), 2) + pos(tri(i,3), 2))/3 ]';
    F = eye(2)+([(e1-E1),(e2-E2)]/[E1,E2])';
   
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %edit the cauchyStress method to add material models
    stress = cauchyStress(F);
    
    %compute forces for each edge
    %edge 1
    n1 = [e1(2) ; -e1(1)];
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set fe1 to the force on the first edge of the triangle (from node 1 to
    %node 2)
    l1 = norm(e1);
    fe1 = stress*n1*l1;
    
    
    f(tri(i,2), :) = f(tri(i,2), :) - 0.5*fe1';
    f(tri(i,1), :) = f(tri(i,1), :) - 0.5*fe1';
    
    %edge2
    n2 = [e3(2) ; -e3(1)];
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set fe2 to the force on the first edge of the triangle (from node 2 to
    %node 3)
    l2 = norm(e3);
    fe2 = stress*n2*l2;
    f(tri(i,3), :) = f(tri(i,3), :) - 0.5*fe2';
    f(tri(i,2), :) = f(tri(i,2), :) - 0.5*fe2';
    
    
    n3 = [-e2(2) ; e2(1)];
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set fe3 to the force on the first edge of the triangle (from node 3 to
    %node 1)
    l3 = norm(e2);
    fe3 = stress*n3*l3;
    f(tri(i,3), :) = f(tri(i,3), :) - 0.5*fe3';
    f(tri(i,1), :) = f(tri(i,1), :) - 0.5*fe3';
    
    
    %distribute mass to all vertices
    m(tri(i,1), : ) = m(tri(i,1), : ) + [me me]/3;
    m(tri(i,2), : ) = m(tri(i,2), : ) + [me me]/3;
    m(tri(i,3), : ) = m(tri(i,3), : ) + [me me]/3;
   
    
end

%reshape force vector and compute accelerations, add gravity here
f = reshape(f', 2*numVerts,1);
m = reshape(m', 2*numVerts,1);
a = f./m + repmat(g, numVerts,1);

end

function status = femOutputFcn(time, state, flag, pos, tri, v)

if strcmp(flag, 'done') == 0
numVerts = size(v,1);
hold on
clf
pos(1:numVerts,:) = [state(1:2:2*numVerts, end) state(2:2:2*numVerts, end)];
posz(1:numVerts) = 0;
trimesh(tri, pos(:,1), pos(:,2));
hold off
drawnow
end

status = 0;
end



function stress = cauchyStress(F)
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set P = to the Piola-Kirchoff 1 stress for a neohookean model.
    %Use parameters mu = 50, lambda = 50 to test
    mu = 50;
    lambda = 50;
%     neo-hookean
    J = det(F);
    P = (mu*(F- inv(F).')) + (lambda*log(J)*inv(F).');
    stress = (P*F')./J;

%     St. Venant-Kirchhoff Model
%     E = 0.5*(F'*F - eye(size(F,1)));
%     P = F*((2*mu*E) + (lambda*trace(E)*eye(size(F,1))));
%     J = det(F);
%     stress = (P*F')./J;

    
end

function [tri, v] = createSquareMesh(width, height, dx, dy)

[X,Y] = meshgrid(1:dx:width, 1:dy:height);

v = [X(:) Y(:)];
tri = delaunay(v(:,1), v(:,2));
disp(tri);
numVerts = size(v,1);
%center mesh
v(:,1) = v(:,1) - sum(v(:,1))./numVerts;
v(:,2) = v(:,2) - min(v(:,2));

end

function [tri,v] = createCircleMesh(r,dtheta)

theta = 0:dtheta:2*pi;
X = r*cos(theta);
Y = r*sin(theta);

v = [X(:) Y(:)];
tri = delaunay(v(:,1) , v(:,2));
end

