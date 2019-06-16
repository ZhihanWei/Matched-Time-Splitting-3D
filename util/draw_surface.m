%********************** Data Clear ***********************%
clear all
clf
format long

%********************** Load File ***********************%
fileID = fopen('Tanglecube_Error_160.txt','r');
%fileID = fopen('Tanglecube_Solution_160.txt','r');
formatSpec = '%f %f %f %12.6f';
A = textscan(fileID,formatSpec);
fclose(fileID);

for i = 1:4
    Info(:,i) = A{i};
end

%********************** Set up mesh ***********************%
nx = 160; dx = 3.99;
ny = 160; dy = 3.99;
nz = 160; dz = 4.99;
X = -dx:(2*dx)/(nx-1):dx;
Y = -dy:(2*dy)/(ny-1):dy;
Z = -dz:(2*dz)/(nz-1):dz;
[x,y,z]=meshgrid(X,Y,Z);

%************** Implicit function definition ************%
w= x.^4-5.*x.^2+y.^4-5.*y.^2+z.^4-5.*z.^2+10;

%********************** Color matrix ***********************%
cdata = ones(size(w))*min(Info(:,4));

%***************** Set up color for surface ******************%
for i = 1:size(Info(:,1))
    for j = 1:4
        cdata(Info(i,1)+1,Info(i,2)+1,Info(i,3)+1) = Info(i,4);
    end 
end

%****************** Draw surface graph *******************%
cdata = smooth3(cdata,'box',7);
p=patch(isosurface(x,y,z,w,0,'noshare'));
cmap = colormap(jet(100));
isonormals(x,y,z,w,p)
isocolors(x,y,z,cdata,p)

%******************* Domain defined *******************%
xlabel('x-axis')
ylabel('y-axis')
zlabel('z-axis')
D = 2.99;
xlim([-D,D])
ylim([-D,D])
zlim([-D,D])

%********************** Custom view ***********************%
p.FaceColor = 'interp';
p.EdgeColor = 'none';
daspect([1 1 1])
camlight
lighting gouraud
view([40,20])
rotate3d on
axis off
zoom(1.6)

%****************** Custom Colorbar *******************%
c = colorbar;  

n = 6;

c.Limits = [min(cdata(:)),max(cdata(:))];
h1 = (max(cdata(:))-min(cdata(:)))/n;
c.Ticks = min(cdata(:)):h1:max(cdata(:));

h2 = (max(Info(:,4))-min(Info(:,4)))/n;
label = min(Info(:,4)):h2:max(Info(:,4));
label = num2str(round(label.*10000)/10000,'%1.4f\n');
set(c,'XTicklabel',label);

c.Location = 'manual';
c.Position = [0.17,0.3,0.03,0.5];
c.Box = 'off';
c.FontSize = 26;
c.Units = 'pixels';
c.FontWeight = 'bold';
c.FontName = 'Arial';
c.TickDirection = 'out';