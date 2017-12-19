% [x,y] = meshgrid(1:15,1:15);
% tri = delaunay(x,y);
% z = peaks(15);
% % z = meshgrid(1:15);
% trimesh(tri,x,y,z)
plot(plot::Sweep([u*cos(u), u*sin(u), u], u = 0..4*PI),
     CameraDirection = [90, 50, 120])