%%
load C-Estimate.mat ; 
% load('E:\Runs\PIG-TWG\C-Estimate.mat')
FC=scatteredInterpolant(xC,yC,C); 
save('FC_jutul_vest.mat','FC')



load AGlen-Estimate.mat ; 
% load('E:\Runs\PIG-TWG\AGlen-Estimate.mat')
FA=scatteredInterpolant(xA,yA,AGlen); 
save('FA_jutul_vest.mat','FA')


