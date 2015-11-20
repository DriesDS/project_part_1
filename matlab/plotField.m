points = load('matlab/plotfieldpoints.out');

x = points(:,1);
y = points(:,2);
z = points(:,3);
figure;
tri = delaunay(x,y);
plot(x,y,'.')

[r,c] = size(tri);
disp(r)

h = trisurf(tri, x, y, z);

print('field', '-dpng')

