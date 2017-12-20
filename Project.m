[tri, v] = import3Dmesh('./obj/untitled.obj');
trimesh(tri,v(:,1), v(:,2), v(:,3));
numVerts = size(v,1);
pos = zeros(numVerts,3);
odeFun = @(time, state)(femOde(time,state,tri,v,pos));
outputFun = @(time, state, flag)(femOutputFcn(time,state,flag,pos, tri, v));
options = odeset('OutputFcn', outputFun);

y0 = [reshape(v', 3*numVerts,1); zeros(3*numVerts,1)];
[time, state] = ode45(odeFun, [0, 10], y0, options);

pos(1:numVerts,:) = [state(end, 1:3:3*numVerts)' state(end, 2:3:3*numVerts)' state(end, 3:3:3*numVerts)'];
trimesh(tri, pos(:,1), pos(:,2), pos(:,3));



function a = femAccelerations(tri, v, pos)

numTris = size(tri, 1);
numVerts = size(v,1);
f = zeros(numVerts,3);
m = zeros(numVerts, 3);
density = 1.0;
g = [0 -9.8 0]';

for i=1:numTris
    
    %Edge vectors
    E1 = (v(tri(i,2), :)-v(tri(i,1), :))';
    E2 = (v(tri(i,3), :)-v(tri(i,1), :))';
    
    e1 = (pos(tri(i,2), :)-pos(tri(i,1), :))';
    e2 = (pos(tri(i,3), :)-pos(tri(i,1), :))';
    e3 = (pos(tri(i,3), :)-pos(tri(i,2), :))';
%     area = 0.5*det([e1 e2]);
    angle = acos(dot(e1,e2)/(norm(e1)*norm(e2)));
    area = 0.5*norm(e1)*norm(e2)*sin(angle);
    me = density*area;
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %compute F
%     disp(pos);
%     cref = [(v(tri(i,1), 1)+v(tri(i,2), 1) + v(tri(i,3), 1))/3 , (v(tri(i,1), 2)+v(tri(i,2), 2) + v(tri(i,3), 2))/3 ]';
%     cpos = [(pos(tri(i,1), 1)+pos(tri(i,2), 1) + pos(tri(i,3), 1))/3 , (pos(tri(i,1), 2)+pos(tri(i,2), 2) + pos(tri(i,3), 2))/3 ]';
    F = eye(3)+([(e1-E1),(e2-E2)]/[E1,E2])';
    disp(E2);
%    F = rand(3);
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %edit the cauchyStress method to add material models
    stress = cauchyStress(F);
    
    %compute forces for each edge
    %edge 1
    n1 = [e1(2) ; e1(1) ; -e1(3)];
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set fe1 to the force on the first edge of the triangle (from node 1 to
    %node 2)
    l1 = norm(e1);
%     disp(e1);
%     disp(n1);
    
    fe1 = stress*n1*l1;
    
    f(tri(i,2), :) = f(tri(i,2), :) - 0.5*fe1';
    f(tri(i,1), :) = f(tri(i,1), :) - 0.5*fe1';
    
    
    %edge2
    n2 = [-e3(2) ; -e3(1) ; e3(3)];
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set fe2 to the force on the first edge of the triangle (from node 2 to
    %node 3)
    l2 = norm(e3);
    fe2 = stress*n2*l2;
    f(tri(i,3), :) = f(tri(i,3), :) - 0.5*fe2';
    f(tri(i,2), :) = f(tri(i,2), :) - 0.5*fe2';
    
    
    n3 = [-e2(2) ; e2(1) ; e2(3)];
    
    %%%%%%%%%%ASSIGNMENT%%%%%%%%%%%%%%
    %Set fe3 to the force on the first edge of the triangle (from node 3 to
    %node 1)
    l3 = norm(e2);
    fe3 = stress*n3*l3;
    f(tri(i,3), :) = f(tri(i,3), :) - 0.5*fe3';
    f(tri(i,1), :) = f(tri(i,1), :) - 0.5*fe3';
    
    
    %distribute mass to all vertices
    m(tri(i,1), : ) = m(tri(i,1), : ) + [me me me]/3;
    m(tri(i,2), : ) = m(tri(i,2), : ) + [me me me]/3;
    m(tri(i,3), : ) = m(tri(i,3), : ) + [me me me]/3;
   
    
end

%reshape force vector and compute accelerations, add gravity here
f = reshape(f', 3*numVerts,1);
m = reshape(m', 3*numVerts,1);
a = f./m + repmat(g, numVerts,1);

end



function status = femOutputFcn(time, state, flag, pos, tri, v)

if strcmp(flag, 'done') == 0
numVerts = size(v,1);
hold on
clf
pos(1:numVerts,:) = [state(1:3:3*numVerts, end) state(2:3:3*numVerts, end) state(3:3:3*numVerts, end)];
trimesh(tri, pos(:,1), pos(:,2), pos(:,3));
hold off
drawnow
end

status = 0;
end

function y = femOde(time, state, tri, v, pos)
    numVerts = size(v,1);
    pos(1:numVerts,:) = [state(1:3:3*numVerts) state(2:3:3*numVerts) state(3:3:3*numVerts)];
    
    %hack to fix object to the floor
    %just set y positions of vertices to zero
    floorVert = abs(v(:,3)) < 1e-6; %some arbitrary toleranc
%     disp(floorVert);
    
    %set positions to zero
    pos(floorVert,3) = 0;
    
    %set velocities to zero
    vel = [state((3*numVerts+1):3:end) state((3*numVerts+2):3:end) state((3*numVerts+3):3:end)];
    vel(floorVert,3) = 0;
    [m,n] = size(state);
%     disp(numVerts);
%     disp(m);
%     disp(n);
    y = [reshape(vel', 3*numVerts,1) ; femAccelerations(tri, v, pos)];
    
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

function [tri,v] = import3Dmesh(filepath)
object = readObj(filepath);
v = [object.v(:,1) object.v(:,2) object.v(:,3)];
tri = object.f;
% disp(teapot.f);
% [m n] = size(tri);
% disp(m);
% disp(n);
% trimesh(tri,v(:,1), v(:,2), v(:,3));
end
