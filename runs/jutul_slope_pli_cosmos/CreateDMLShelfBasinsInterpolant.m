basins = s.* 0;
mask_ekstrom = MUA.coordinates(:,1) > -400e3 & MUA.coordinates(:,1) <= -220e3;
mask_jelbart = MUA.coordinates(:,1) > -220e3 & MUA.coordinates(:,1) <= -100e3;
mask_jutulst = MUA.coordinates(:,1) > -100e3 & MUA.coordinates(:,1) <=  200e3;

basins(mask_ekstrom) = 1;
basins(mask_jelbart) = 2;
basins(mask_jutulst) = 3;

figure; 
PlotMeshScalarVariable(CtrlVar,MUA,basins); 
hold on; 
[xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');

DMLBasinsInterpolant=scatteredInterpolant(MUA.coordinates(:,1),MUA.coordinates(:,2),basins,'nearest','nearest');
save('DMLBasinsInterpolant')