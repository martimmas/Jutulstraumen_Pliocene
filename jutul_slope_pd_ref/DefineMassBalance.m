function [UserVar,as,ab,dasdh,dabdh]=DefineMassBalance(UserVar,CtrlVar,MUA,time,s,b,h,S,B,rho,rhow,GF)

persistent Fas Fab Fas_lgm Ft_clim Fp_clim Fsr

if contains(UserVar.RunType,'Forward') ==1 % i.e., diagnostic and transient
    f=1;
    gammaT = 1e-4;
    rhow = 1028; 
    cp = 3974;
    rhoi = 917;
    Li = 3.34e5;

    % t0 and s0 values for PICO should be averaged along the ice front at z = b
    s0 = 34.48; % From Reese et al. (2018) - Schmidtko et al. (2014)
    if ~isfield(UserVar.Qbm,'Temp0')
        t0 = -1.65 ; % -1.65 From Reese et al. (2018) - Schmidtko et al. (2014)
    else
        t0 = UserVar.Qbm.Temp0;
    end

    l1=-0.0575; l2=0.0832; l3=7.59e-4;
    tf = l1*s0 + l2 + l3*(b-S);
    %tf(h <= CtrlVar.ThickMin & GF.node <1) = NaN; % CANNOT ADD NANS TO IT


    if UserVar.SMB.ReadFromFile && isempty(Fas)
        fprintf('SurfaceMassBalance: loading file: %-s ',UserVar.SMB.File)
        Fas = load(UserVar.SMB.File).F;
        fprintf(' done \n')
        %as=Fas(MUA.coordinates(:,1),MUA.coordinates(:,2));
        %as(abs(as)>1e10) = 0;
        if isfield(UserVar.SMB,'TargetFile') ==1
            Fas_lgm = load(UserVar.SMB.TargetFile).F;
        end
    end
    
    if UserVar.SMB.ComputeFromTS && (isempty(Ft_clim) || isempty(Fp_clim) || isempty(Fsr))
        Ft_clim = load(UserVar.SMB.Tclim_file).Ft_clim;
        Fp_clim = load(UserVar.SMB.Pclim_file).Fp_clim;
        Fsr     = load(UserVar.SMB.s_ref_file).Fsr;
    end

    if strcmp(UserVar.Qbm.Mode,'ReadFromFile')==1 && isempty(Fab)
        fprintf('SurfaceMassBalance: loading file: %-s ',UserVar.Qbm.File)
        Fab = load(UserVar.Qbm.File).F;
        fprintf(' done \n')
        %ab=Fab(MUA.coordinates(:,1),MUA.coordinates(:,2));
        %ab(abs(ab)>1e10) = 0;
    end
        
    if (UserVar.SMB.ReadFromFile ~= 1) 
        as=0 ;
    else
        as=Fas(MUA.coordinates(:,1),MUA.coordinates(:,2));                        
        
        if ((isfield(UserVar.SMB,'TimeEnd')==1 && (time >= UserVar.SMB.TimeInit))) && ...
            ((isfield(UserVar.SMB,'TargetFile')==1) || UserVar.SMB.ComputeFromTS==1 )...
            
            
            if UserVar.SMB.ComputeFromTS
        
                % create climatology matrices and populate them with interpolants
                temp_clim = zeros(length(MUA.coordinates(:,1)),12); 
                prec_clim = zeros(length(MUA.coordinates(:,1)),12);
                for n=1:12
                    temp_clim(:,n) = Ft_clim{n}(MUA.coordinates(:,1),MUA.coordinates(:,2));
                    prec_clim(:,n) = Fp_clim{n}(MUA.coordinates(:,1),MUA.coordinates(:,2));
                end
                s_ref = Fsr(MUA.coordinates(:,1),MUA.coordinates(:,2)); % reference elevation for LGM outputs
                as_lgm = Calc_SMB_from_atm_flds(UserVar,s,temp_clim,prec_clim,s_ref);
            else
                as_lgm = Fas_lgm(MUA.coordinates(:,1),MUA.coordinates(:,2));
            end            
            
            as_anom_rate = (as_lgm - as)/(UserVar.SMB.TimeEnd-UserVar.SMB.TimeInit);
            as = as + as_anom_rate*(time-UserVar.SMB.TimeInit);
            
            if time >= UserVar.SMB.TimeEnd
                as = as_lgm;
            end
        end
        
        as(abs(as)>1e4) = 0;    
    end

    if strcmp(UserVar.Qbm.Mode,'ReadFromFile') ~= 1
        ab=0;    
    end

    switch UserVar.Qbm.Mode

        case 'ReadFromFile'
            ab=-1*Fab(MUA.coordinates(:,1),MUA.coordinates(:,2));
            ab(abs(ab)>1e10) = 0;        
            ab=ab.*(1-GF.node); % applies only to floating nodes

        case 'Quadratic'                        
            M = f .* gammaT .* ((rhow * cp)./(rhoi * Li)).^2 * (t0 - tf) .* abs(t0-tf);
            ab = -M * 365.25 * 86400;
            ab=ab.*(1-GF.node); % applies only to floating nodes

        case 'PICO'
            PICO_opts.InfoLevel=0;            
            %PICO_opts.MeshBoundaryCoordinates = CreateMeshBoundaryCoordinatesForJutulstraumen(CtrlVar,'pd');
            %load('DMLBasinsInterpolant','DMLBasinsInterpolant');
            %PICO_opts.BasinsInterpolant = DMLBasinsInterpolant;
            
            %PICO_opts.ContinentArea = 1e9; % default: 5e10;
            PICO_opts.nmax = 6;
            PICO_opts.C1 = 1e6;
            PICO_opts.gamTstar = 7.4074e-5; % This should be the maximum allowed!! default was 2e-5
            PICO_opts.algorithm = 'watershed'; % test 'watershed', 'polygon' or 'oneshelf'             
            PICO_opts.Tbasins = t0; %[t0-0.1,t0+0.2,t0];
            PICO_opts.Sbasins = s0; %[s0, s0, s0];
            ab = PICO(UserVar,CtrlVar,MUA,GF,h,rhoi,rhow,PICO_opts);
            
            if UserVar.Qbm.GLenhance==1            
                cooA=[MUA.coordinates(:,1) MUA.coordinates(:,2)];
                KdTree=KDTreeSearcher(cooA);
                CtrlVar.PlotGLs=0; CtrlVar.GLsubdivide=1; CtrlVar.LineUpGLs=0;
                [xGL,yGL]=PlotGroundingLines(CtrlVar,MUA,GF); % no need to align GL.
                ds = 5000;
                [ID,~,~,~]=FindAllNodesWithinGivenRangeFromGroundingLine(CtrlVar,MUA,xGL,yGL,ds,KdTree);
                ab(ID) = ab(ID) .* UserVar.Qbm.GLenhanceFactor;
            end
            
        case 'none' 
            ab= 0;
    end
    
    % Taking care that are no NaN's in the field
    if any(isnan(ab))
        ab(isnan(ab))=0.0; as(isnan(as))=0.0;
    end
    
else
    as=0;
    ab=0;
end
dasdh=0;
dabdh=0;

if isfield(UserVar.Qbm,'Calving')
    if UserVar.Qbm.Calving==1
    
    % Here a fictitious basal melt distribution is applied over the ice shelf downstream
    % of x=400km for the first few years to melt away all/most floating ice. 
    %
    % The melt is prescribed as a function of ice thickness and to speed things up
    % the mass-balance feedback is provided here as well. This requires setting 
    %
    %   CtrlVar.MassBalanceGeometryFeedback=3;
    %
    % in DefineInitialInputs.m
    %
    
    dabdh=zeros(MUA.Nnodes,1) ;
    dasdh=zeros(MUA.Nnodes,1) ;
    
%     if (CtrlVar.time+CtrlVar.dt)  < 5
        
        GF=IceSheetIceShelves(CtrlVar,MUA,GF);
        NodesSelected=(h<UserVar.CalvingThreshold) & GF.NodesDownstreamOfGroundingLines;
        
        % ab = -(h-hmin)  , dab=-1 ;
        ab(NodesSelected)=-(h(NodesSelected)-CtrlVar.ThickMin) ;
        dabdh(NodesSelected)=-1;
        
        
%     end
    end
end

end