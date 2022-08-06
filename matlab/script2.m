%%
nb=10000;
t=1;
da=1.0;
db=0.5;%0.5;
k=0.059;
f=0.042;
[V,F]=readOBJ('../growing-shapes/spot.obj');
o=V(100,:);
for i=1:2
    [V,F]=loop(V,F);
end
L=cotmatrix(V,F);
% M=massmatrix(V,F);
% L=M\L;
a=ones(size(L,1),1);
b=zeros(size(L,1),1);
for i=1:size(V,1)
    if(norm(o-V(i,:)))<0.01
        b(i)=1;
    end
end
%%
%the size of the motif depends on how smooth the mesh is (more triangles with smaller areas, more motifs)
first=speye(size(L,1))-t*da*L;
%first=first*first*first*first;
second = speye(size(L,1))-t*db*L;
%second=second*second*second*second;
for i=1:nb
    i
    termbyterm = a.*b.*b;
    a = first\(a + t*(-termbyterm+f*(1-a)));
    b = second\(b + t*(termbyterm-(k+f)*b));
    grid=(a)*255;
    grid=min(max(grid,0),255);
    %image(grid);
    tsurf(F,V,'CData',grid);shading interp;
    pause(1);
end