
function MeshBoundaryCoordinates=CreateMeshBoundaryCoordinates(CtrlVar,fin)


% load PigBoundaryCoordinates MeshBoundaryCoordinates
% MeshBoundaryCoordinates=flipud(MeshBoundaryCoordinates);
% I=MeshBoundaryCoordinates(:,2)>-973818.250320996-1;
% MeshBoundaryCoordinates=[MeshBoundaryCoordinates(I,1) MeshBoundaryCoordinates(I,2)];

% grd_coords = load('/datadisk/phd/Jutul/mesh_creation/bc_grd_jutulstraumen.dat');
grd_coords = load(fin);


mnt_coords = grd_coords;

x0 = mnt_coords(:,1); y0 = mnt_coords(:,2);


CtrlVar.GLtension=1e-9; % tension of spline, 1: no smoothing; 0: straight line
CtrlVar.GLds=CtrlVar.MeshSizeMin;  %CtrlVar.MeshSizeBoundary;

[xB,yB,nx,ny] = Smooth2dPos(x0,y0,CtrlVar);
CtrlVar.GLtension=1;
[xB,yB,nx,ny] = Smooth2dPos(xB,yB,CtrlVar);

% xB=x0;
% yB=y0;

%figure; plot(xB,yB)
MeshBoundaryCoordinates=[xB,yB];

%CtrlVar.GLtension=1e-9; % tension of spline, 1: no smoothing; 0: straight line
%CtrlVar.GLds=CtrlVar.MeshSizeBoundary;
%[x,y,nx,ny] = Smooth2dPos(jutul_coords(:,1),jutul_coords(:,2),CtrlVar);
%MeshBoundaryCoordinates=[x(:) y(:)];


end