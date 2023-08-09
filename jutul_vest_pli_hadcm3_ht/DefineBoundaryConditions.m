function  [UserVar,BCs]=DefineBoundaryConditions(UserVar,CtrlVar,MUA,BCs,time,s,b,h,S,B,ub,vb,ud,vd,GF)


persistent AA BB %Af Bf


if isempty(AA)
  
    % load points that define the line-segments along which the BCs are to be defined
    grd_coords = load([UserVar.IOFilesPath,'/MeshFiles/mesh_lines_jutul_vest/bc_grd_vest_jutul.dat']);
    xx = grd_coords(:,1); yy = grd_coords(:,2);
    AA=[xx(1:end-1) yy(1:end-1)] ; BB=[xx(2:end) yy(2:end)];

end


% find all boundary nodes within 1m distance from the line segment.
x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2); tolerance=CtrlVar.MeshSizeMax*2;
I = DistanceToLineSegment([x(MUA.Boundary.Nodes) y(MUA.Boundary.Nodes)],AA,BB,tolerance);
BCs.vbFixedNode=MUA.Boundary.Nodes(I);
BCs.ubFixedNode=MUA.Boundary.Nodes(I);

% Implement zero velocities at the inland domain boundaries
BCs.ubFixedValue=BCs.ubFixedNode*0;
BCs.vbFixedValue=BCs.vbFixedNode*0;

% Implement measures velocities at the inland domain boundaries
% [BCs.ubFixedValue,BCs.vbFixedValue]=EricVelocities(CtrlVar,[x(MUA.Boundary.Nodes(I)) y(MUA.Boundary.Nodes(I))]);

%% Implement measures velocities at the floating end
% if isempty(Af)
%     flt_coords = load('mesh_lines/bc_flt_jutulstraumen.dat');
%     xx = flt_coords(:,1); yy = flt_coords(:,2);
%     Af=[xx(1:end-1) yy(1:end-1)] ; Bf=[xx(2:end) yy(2:end)];
% end
% I = DistanceToLineSegment([x(MUA.Boundary.Nodes) y(MUA.Boundary.Nodes)],Af,Bf,tolerance);
% ubFltNodes=MUA.Boundary.Nodes(I);
% vbFltNodes=MUA.Boundary.Nodes(I);
% [ubFltValue,vbFltValue]=EricVelocities(CtrlVar,[x(ubFltNodes) y(vbFltNodes)]);
% 
% 
% BCs.vbFixedNode=[BCs.vbFixedNode; vbFltNodes];
% BCs.ubFixedNode=[BCs.ubFixedNode; ubFltNodes];
% 
% BCs.ubFixedValue=[BCs.ubFixedValue; ubFltValue];
% BCs.vbFixedValue=[BCs.vbFixedValue; vbFltValue];

end