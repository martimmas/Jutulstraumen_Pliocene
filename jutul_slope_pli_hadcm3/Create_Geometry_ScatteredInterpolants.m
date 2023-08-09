% Creates a GeometryInterpolant from a Ua RestartFile
run = 'Jutul_picorac_ctrl_rlx9-Forward-Transient-0000100';

load(['ResultsFiles/',run])
x = MUA.coordinates(:,1); y = MUA.coordinates(:,2);

Fs=scatteredInterpolant(x,y,s); 
FB=scatteredInterpolant(x,y,B); 
Fh=scatteredInterpolant(x,y,h); 
save(['ResultsFiles/GeometryInterpolant_',run,'.mat'],'Fs','FB','Fh')

