function [UserVar,AGlen,n]=DefineAGlenDistribution(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)


persistent FA

n=3;

if ~UserVar.AGlen.ReadFromFile    
         
    if UserVar.Temp.ReadFromFile
        
        if isempty(FA)
            load('InputAC/tempavg_2d_a08_pd-racmo_calib-sld-ocn_all-vars_002kyr_Temp2D.mat','F');
            FA = F;
        end
        
        temp = FA(MUA.coordinates(:,1),MUA.coordinates(:,2));
        AGlen=AGlenVersusTemp(temp);
    else
        AGlen=AGlenVersusTemp(-10);
    end
    
else
    
    if isempty(FA)
        
        if isfile(UserVar.AFile)
            
            fprintf('DefineSlipperyDistribution: loading file: %-s ',UserVar.AFile)
            load(UserVar.AFile,'FA')
            fprintf(' done \n')
            
        else
            
            load('AGlen-Estimate.mat','AGlen','xA','yA')
            FA=scatteredInterpolant(xA,yA,AGlen);
            save(UserVar.AFile,'FA')
            
        end
        
    end
    
    AGlen=FA(MUA.coordinates(:,1),MUA.coordinates(:,2));        
    
end

