%%
load('C-Estimate.mat')
FC=scatteredInterpolant(xC,yC,C); 
save('InputAC/FC_inv_sico_uvdhdt_err50cm.mat','FC')

load('AGlen-Estimate.mat')
FA=scatteredInterpolant(xA,yA,AGlen); 
save('InputAC/FA_inv_sico_uvdhdt_err50cm.mat','FA')


