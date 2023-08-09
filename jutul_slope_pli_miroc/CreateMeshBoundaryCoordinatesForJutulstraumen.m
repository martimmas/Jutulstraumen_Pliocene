
function MeshBoundaryCoordinates=CreateMeshBoundaryCoordinatesForJutulstraumen(CtrlVar,ext)


% load PigBoundaryCoordinates MeshBoundaryCoordinates
% MeshBoundaryCoordinates=flipud(MeshBoundaryCoordinates);
% I=MeshBoundaryCoordinates(:,2)>-973818.250320996-1;
% MeshBoundaryCoordinates=[MeshBoundaryCoordinates(I,1) MeshBoundaryCoordinates(I,2)];

% grd_coords = load('/datadisk/phd/Jutul/mesh_creation/bc_grd_jutulstraumen.dat');
grd_coords = load('mesh_lines/bc_grd_jutulstraumen.dat');
if strcmp(ext,'extended') == 1
    flt_coords = load('mesh_lines/bc_flt_ext_jutulstraumen.dat'); % encompasses LGM maximum extent
else
    flt_coords = load('mesh_lines/bc_flt_jutulstraumen.dat');    
end

jutul_coords = [grd_coords; flt_coords];

x0 = jutul_coords(:,1); y0 = jutul_coords(:,2);

CtrlVar.GLtension=1e-9; % tension of spline, 1: no smoothing; 0: straight line
CtrlVar.GLds=CtrlVar.MeshSizeMin; %CtrlVar.MeshSizeBoundary;

[xB,yB,nx,ny] = Smooth2dPos(x0,y0,CtrlVar);

%xB=x0;
%yB=y0;

MeshBoundaryCoordinates=[xB,yB];

%CtrlVar.GLtension=1e-9; % tension of spline, 1: no smoothing; 0: straight line
%CtrlVar.GLds=CtrlVar.MeshSizeBoundary;
%[x,y,nx,ny] = Smooth2dPos(jutul_coords(:,1),jutul_coords(:,2),CtrlVar);
%MeshBoundaryCoordinates=[x(:) y(:)];


end