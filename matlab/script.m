%%
nb=10000;
t=1;
da=1.0;
db=0.5;%0.5;
k=0.06;
f=0.025;
grid_a = ones(80);
grid_b = zeros(80);
% for i=1:20
%     grid_b(x(i),y(i))=1;
% end
% grid_b(400:520,500:520)=1;
% grid_b(575:695,500:520)=1;
%l=100:80:900;
grid_b(40:60,40:60)=1;
%grid_a(40:60,40:60)=0;
L = [0.05 0.2 0.05;0.2 -1 0.2;0.05 0.2 0.05];
%L=[0 0.25 0;0.25 -1 0.25;0 0.25 0];
%%

for i=1:nb
    i
    c1 = conv2(grid_a,L,'same');
    c2 = conv2(grid_b,L,'same');
    termbyterm = grid_a.*grid_b.*grid_b;
    new_grid_a = grid_a + t*(da*c1-termbyterm+f*(1-grid_a));
    new_grid_b = grid_b + t*(db*c2+termbyterm-(k+f)*grid_b);
    grid_a = new_grid_a;
    grid_b = new_grid_b;
    grid=(grid_a)*255;
    grid=min(max(grid,0),255);
    image(grid);
    pause(0.1);
end