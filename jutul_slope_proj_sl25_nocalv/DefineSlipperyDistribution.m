function [UserVar,C,m,q,muk]=DefineSlipperyDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

persistent FC

m=3;

if ~UserVar.Slipperiness.ReadFromFile
            
    ub=10 ; %originally 10
    tau=80 ; % units meters, year , kPa
    C0=ub/tau^m;
    C=C0;
    
else
    
    
    if isempty(FC)
        
        if isfile(UserVar.CFile)
            fprintf('DefineSlipperyDistribution: loading file: %-s ',UserVar.CFile)
            load(UserVar.CFile,'FC')
            fprintf(' done \n')
        else
            % create a FC file
            load('C-Estimate.mat','C','xC','yC')
            FC=scatteredInterpolant(xC,yC,C); 
            save(UserVar.CFile,'FC')
        end
    end
    
    C=FC(MUA.coordinates(:,1),MUA.coordinates(:,2));        
    
end

% m=3; % defined up top!
q=1 ;      % only needed for Budd sliding law
muk=0.5 ;  % required for Coulomb friction type sliding law as well as Budd, minCW (Tsai), rCW  (Umbi) and rpCW (Cornford).    
end