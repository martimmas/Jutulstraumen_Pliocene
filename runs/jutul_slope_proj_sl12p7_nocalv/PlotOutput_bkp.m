function  PlotOutput(InputFile,plots)
addpath /mnt/3tb/Ua/cmocean
addpath /mnt/3tb/Ua/cbrewer
addpath ../../

load(InputFile);

if exist('F','var') == 1
    v2struct(struct(F));
end

if exist('CtrlVarInRestartFile','var') == 1
    CtrlVar = CtrlVarInRestartFile;
end

% if exist('UserVarInRestartFile','var') == 1
%     UserVar = UserVarInRestartFile;
% end


if strcmp(plots,'')
    CtrlVar.DefineOutputs='-sbB-ubvb-s-save-';
else
    CtrlVar.DefineOutputs = plots;
end

%%
% if ~isfield(CtrlVar,'DefineOutputs')
%     CtrlVar.uvPlotScale=[];
%     %plots='-ubvb-udvd-log10(C)-log10(Surfspeed)-log10(DeformationalSpeed)-log10(BasalSpeed)-log10(AGlen)-';
%     plots='-ubvb-log10(BasalSpeed)-sbB-ab-log10(C)-log10(AGlen)-';    
% else
%     plots=CtrlVar.DefineOutputs;
% end

CtrlVar.QuiverColorPowRange=2; 

%%
GLgeo=GLgeometry(MUA.connectivity,MUA.coordinates,GF,CtrlVar);
TRI=[]; DT=[]; xGL=[] ; yGL=[] ;
x=MUA.coordinates(:,1);  y=MUA.coordinates(:,2);
%% Plots

if contains(plots,'-B')
    figure('Renderer','painters','Position',[10 800 600 500]);
    PlotMeshScalarVariable(CtrlVar,MUA,B);
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')    
    xlabel('x (km)') ; ylabel('y (km)')
    cbh = colorbar; title(cbh,'(m)')
    zlimits = [min(B) max(B)];
    cmapsea = flip(cbrewer('seq','Blues',32));
    cmapsea = cmapsea(1:20,:);
    demcmap(zlimits,64,cmapsea,[]); % remove white here? how?
    box on;
    if contains(plots,'-Bgl-')
        xlim([-370 125])
        ylim([1850 2250])    
    end
    if contains(plots,'-savef-')
        f = gcf;                
        set(f,'color','w')
        expname = strsplit(InputFile,'/'); expname = expname{end};
        fig_name = ['Figures/Bedrock-',expname];
        export_fig(fig_name,'-png','-preserve_size','-r300')
    else
        title(sprintf('Bedrock at t=%-g ',CtrlVar.time)) ;
    end
end

%%
if contains(plots,'-s-')
    splot = s;
    %splot(splot-b<=CtrlVar.ThickMin) = NaN;
    figure('Renderer','painters','Position',[10 800 600 500]);
    set(gca,'FontSize',12)
    PlotMeshScalarVariable(CtrlVar,MUA,splot);
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')    
    cbh = colorbar; title(cbh,'(m a.s.l.)')
    xlabel('x (km)') ; ylabel('y (km)')    
    caxis([0 3500])    
    colormap(cmocean('-ice'))
    box on;
    if contains(plots,'-savef-')
        f = gcf;                
        set(f,'color','w')
        expname = strsplit(InputFile,'/'); expname = expname{end};
        fig_name = ['Figures/Surface-',expname];
        export_fig(fig_name,'-png','-preserve_size','-r300')
    else
        title(sprintf('Ice surface at t=%-g ',CtrlVar.time)) ;
    end
end

%%
if contains(plots,'-SurfSpeed-')    
    us=ub+ud;  vs=vb+vd;
    SurfSpeed=sqrt(us.*us+vs.*vs);    
    SurfSpeed(SurfSpeed<=0) = 1e-4;
    figure('Renderer','painters','Position',[10 800 600 500]);
    set(gca,'fontsize',12)
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,SurfSpeed,CtrlVar);
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k','linewidth',1.4);
    PlotMuaBoundary(CtrlVar,MUA,'k')
%     caxis([0 3])    
    [cmap]=cbrewer('div', 'Spectral', 255, 'linear');
    colormap(flip(cmap))
     xlabel('x (km)') ; ylabel('y (km)')
    cbh = colorbar;
    title(cbh,'m/yr')
%     cb_ticks = cbh.Ticks;
%     vel_ticks = 10 .^cb_ticks;
%     cbh.TickLabels=compose('%d',int64(vel_ticks));
    set(gca,'ColorScale','log')
    caxis([0.1 1000])
    box on

    if contains(plots,'-reg-')
        I=h<=10;%CtrlVar2.ThickMin;
        scatter(MUA.coordinates(I,1)/CtrlVar.PlotXYscale,MUA.coordinates(I,2)/CtrlVar.PlotXYscale,100,'.m')
        ylim([1750 2000]); xlim([-150 150])
    end

    if contains(plots,'-savef-')
        f = gcf;                
        set(f,'color','w')
        expname = strsplit(InputFile,'/'); expname = expname{end};
        fig_name = ['Figures/SurfSpeed-',expname];
        export_fig(fig_name,'-png','-preserve_size','-r300')
    else
        title(sprintf('Surface speed at t=%-g ',CtrlVar.time)) ;
    end
end

%%
if contains(plots,'-ab-')
    figure('Renderer','painters','Position',[10 800 600 500])
    set(gca,'FontSize',12)        
    PlotMeshScalarVariable(CtrlVar,MUA,-ab)
    hold on
    %PlotMeshScalarVariable(CtrlVar,MUA,ab_mlt)
    %plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')    
    cbh = colorbar;
    title(cbh,'m/yr')
    xlabel('x (km)') ; ylabel('y (km)') ;% title(colorbar,'(m/yr)')            
    colormap(redblueTecplot);
    caxis([-10 10])    
    %cbh.TickLabels=compose('%d',int64([0., 0.1, 1.0, 2.0]));    
    xlim([-370 125])
    ylim([1850 2250])
    box on
    %axis equal
    if contains(plots,'-savef-')
        f = gcf;                
        set(f,'color','w')
        expname = strsplit(InputFile,'/'); expname = expname{end};
        fig_name = ['Figures/Qbm-',expname];
        export_fig(fig_name,'-png','-preserve_size','-r300')
    else
        title(sprintf('Basal melting at t=%-g ',CtrlVar.time)) ; 
    end
end
%%

if contains(plots,'-logab-')
    
    % Plot the freezing values
    figure('Renderer','painters','Position',[10 800 700 520])
    set(gca,'FontSize',12)    
    ab_freeze = ab;
    ab_freeze(ab <=0) = nan;
    [fh,cbh] = PlotMeshScalarVariable(CtrlVar,MUA,ab_freeze);
    hold on
    %PlotMeshScalarVariable(CtrlVar,MUA,ab_mlt)
    %plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')    
    %cbh = colorbar('east'); 
    set(cbh,'position',[.92 .125 .03 .4],'fontsize',10)    
    set(gca,'ColorScale','log')
    xlabel('x (km)') ; ylabel('y (km)') ;% title(colorbar,'(m/yr)')        
    cmp = redblueTecplot;
    colormap(flip(cmp(1:end/2,:)));
    caxis([0 3])    
    cbh.Ticks=[0., 0.5, 1.0, 2.0];    
    cbh.TickLabels=compose('-%0.1f',[0., 0.5, 1.0, 2.0]);    
    set(cbh, 'YDir', 'reverse' );
    xlim([-370 125])
    ylim([1850 2250])
    box on
    %axis equal
    if contains(plots,'-savef-')
        f = gcf;                
        set(f,'color','w')
        expname = strsplit(InputFile,'/'); expname = expname{end};
        fig_name = ['Figures/Qbm_freeze-',expname];
        export_fig(fig_name,'-png','-preserve_size','-transparent','-r300')
        
    else
        title(sprintf('Basal freezing at t=%-g ',CtrlVar.time)) ; 
    end
    
     % Plot the melting values
    figure('Renderer','painters','Position',[10 800 700 520])    
    set(gca,'FontSize',12)
    ab_melt = ab;
    ab_melt(ab >=0) = nan;
    [fhm,cbhm] = PlotMeshScalarVariable(CtrlVar,MUA,-ab_melt);
    hold on    
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')            
    set(cbhm,'position',[.92 .525 .03 .4],'fontsize',10)
    title(cbhm,'m/yr')
    set(gca,'ColorScale','log')
    xlabel('x (km)') ; ylabel('y (km)') ;% title(colorbar,'(m/yr)')        
    cmp = flip(redblueTecplot);
    colormap(flip(cmp(1:end/2,:)));
    caxis([0 25])    
    cbhm.Ticks=[0., 0.5, 1.0, 2.0]*10;    
    cbhm.TickLabels=compose('%d',int64([0., 0.5, 1.0, 2.0]*10));    
%     set(cbh, 'YDir', 'reverse' );
    xlim([-370 125])
    ylim([1850 2250])
    box on
    %axis equal
    if contains(plots,'-savef-')
        f = gcf;                
        set(f,'color','w')
        expname = strsplit(InputFile,'/'); expname = expname{end};
        fig_name = ['Figures/Qbm_melt-',expname];
        export_fig(fig_name,'-png','-preserve_size','-transparent','-r300')
    else
        title(sprintf('Basal melting at t=%-g ',CtrlVar.time)) ; 
    end
end

%%
if contains(plots,'-as-')

    figure('Renderer','painters','Position',[10 800 600 500])
    set(gca,'FontSize',12)
    PlotMeshScalarVariable(CtrlVar,MUA,as)    
    hold on
    %plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')
    xlabel('x (km)') ; ylabel('y (km)') ; title(colorbar,'(m/yr)')
    colormap(flip(redblueTecplot));
    caxis([-1. 1.])
    axis equal
    box on
    if contains(plots,'-savef-')
        f = gcf;                
        set(f,'color','w')
        expname = strsplit(InputFile,'/'); expname = expname{end};
        fig_name = ['Figures/SMB-',expname];
        export_fig(fig_name,'-png','-preserve_size','-r300')
    else
        title(sprintf('SMB at t=%-g ',CtrlVar.time)) ; 
    end
end

%%
if contains(plots,'-h-')

    figure('Renderer','painters','Position',[10 800 600 500])
    set(gca,'FontSize',12)
    PlotMeshScalarVariable(CtrlVar,MUA,h)
    hold on
    plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
    
    I=h<=10;
    plot(MUA.coordinates(I,1)/CtrlVar.PlotXYscale,MUA.coordinates(I,2)/CtrlVar.PlotXYscale,'.r')
    xlabel('x (km)') ; ylabel('y (km)') ; title(colorbar,'(m/yr)')
    axis equal
    box on;
    if contains(plots,'-savef-')
        f = gcf;          
        set(f,'color','w')
        expname = strsplit(InputFile,'/'); expname = expname{end};
        fig_name = ['Figures/Thickness-',expname];
        export_fig(fig_name,'-png','-preserve_size','-r300')
    else
        title(sprintf('Ice thickness at t=%-g ',CtrlVar.time)) ;
    end
end


%%=================================================================================================
% Plots that might be requested later
%%=================================================================================================

%%    
if contains(plots,'-log10(AGlen)-')
    figure
    PlotMeshScalarVariable(CtrlVar,MUA,log10(AGlen));
    hold on ; 
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')
    title(sprintf('log_{10}(AGlen) at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')
    title(colorbar,'log_{10}(yr^{-1} kPa^{-3})')
    caxis([-9 -7.5])
    box on;
end

%%    
if contains(plots,'-log10(C)-')
    figure
    PlotMeshScalarVariable(CtrlVar,MUA,log10(C));
    hold on 
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')
    title(sprintf('log_{10}(C) at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')
    title(colorbar,'log_{10}(m yr^{-1} kPa^{-3})')
    caxis([-7 1])
    box on;
end


if contains(plots,'-m-')
    figure
    PlotMeshScalarVariable(CtrlVar,MUA,m);
    hold on 
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')
    title(sprintf('Sliding law m at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')
    title(colorbar,'m')    
    box on;
end

if contains(plots,'-n-')
    figure
    PlotMeshScalarVariable(CtrlVar,MUA,n);
    hold on 
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')
    title(sprintf('Glens law n at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')
    title(colorbar,'m')    
    box on;
end


%%
if contains(plots,'-b-')    
    figure ;
    PlotMeshScalarVariable(CtrlVar,MUA,b);
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'k')
    title(sprintf('Ice base at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')
    cbh = colorbar; title(cbh,'(m)')
end





%%=================================================================================================
% Other auxiliary plots
%%=================================================================================================
%%
if contains(plots,'-mesh-')
    figure;
    CtrlVar.MeshColor='b';
    CtrlVar.NodeColor='b';
    PlotMuaMesh(CtrlVar,MUA);
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'color',[0.9290, 0.6940, 0.1250]);
    PlotMuaBoundary(CtrlVar,MUA,'k')    
end

%%
if contains(plots,'-BCs-')    
    figure ;
    PlotBoundaryConditions(CtrlVar,MUA,BCs)
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'r');
    PlotMuaBoundary(CtrlVar,MUA,'b')    
end


%%
if contains(plots,'-ubvb-')
    % plotting horizontal velocities
    figure
    N=1;
    %speed=sqrt(ub.*ub+vb.*vb);
    %CtrlVar.MinSpeedWhenPlottingVelArrows=0; CtrlVar.MaxPlottedSpeed=max(speed); %
    CtrlVar.VelPlotIntervalSpacing='log10';
    %CtrlVar.VelColorMap='hot';
    %CtrlVar.RelativeVelArrowSize=10;
    
    QuiverColorGHG(x(1:N:end),y(1:N:end),ub(1:N:end),vb(1:N:end),CtrlVar);
    hold on ; 
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'b')
    title(sprintf('(ub,vb) at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')  
end


%%
if contains(plots,'-udvd-')
    % plotting horizontal velocities
    figure
    N=1;
    %speed=sqrt(ud.*ud+vd.*vd);
    %CtrlVar.MinSpeedWhenPlottingVelArrows=0; CtrlVar.MaxPlottedSpeed=max(speed); 
    CtrlVar.VelPlotIntervalSpacing='log10';
    %CtrlVar.RelativeVelArrowSize=10;
    %CtrlVar.VelColorMap='hot';
    QuiverColorGHG(x(1:N:end),y(1:N:end),ud(1:N:end),vd(1:N:end),CtrlVar);
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'b')
    title(sprintf('(ud,vd) at t=%-g ',CtrlVar.time)) ; xlabel('xps (km)') ; ylabel('yps (km)')
    axis equal tight
    
end


%%
if contains(plots,'-ub-')  
    figure
    %[FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,ub,CtrlVar)    ;
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,ub,CtrlVar)    ;
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'b')
    title(sprintf('ub at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')
end



%%
if contains(plots,'-log10(SurfSpeed)-')
    
    us=ub+ud;  vs=vb+vd;
    SurfSpeed=sqrt(us.*us+vs.*vs);
    
    figure
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(SurfSpeed),CtrlVar);
    hold on ;
    [xGL,yGL,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'b')        
    title(sprintf('log_{10}(Surface speed) at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')
    cbh = colorbar;
    title(cbh,'log_{10}(m/yr)')    
end


%%
if contains(plots,'-log10(BasalSpeed)-')
    BasalSpeed=sqrt(ub.*ub+vb.*vb); 
    figure
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(BasalSpeed),CtrlVar);
    hold on
    plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
    title(sprintf('log_{10}(Basal speed) at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)') ; title(colorbar,'log_{10}(m/yr)')
end


%%
if contains(plots,'-log10(DeformationalSpeed)-')
    DeformationalSpeed=sqrt(ud.*ud+vd.*vd); 
    figure
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,log10(DeformationalSpeed),CtrlVar);
    hold on
    plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
    title(sprintf('log_{10}(Deformational speed) at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)') ; title(colorbar,'log_{10}(m/yr)')    
end


%%
if contains(plots,'-C-')
    figure
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,C,CtrlVar);
    title(sprintf('C at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')    
    set(gca,'ColorScale','log')
end

if contains(plots,'-AGlen-')
    figure
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,AGlen,CtrlVar);
    title(sprintf('AGlen at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')    
    set(gca,'ColorScale','log')
end

if contains(plots,'-dhdt-')
    figure
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,dhdt,CtrlVar);
    title(sprintf('dh/dt at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')   
    bnds = 0.75*nanmax(abs(dhdt)); caxis([-bnds bnds]);
    colormap(redblueTecplot);    
end

%%
if contains(plots,'-e-')
    % plotting effectiv strain rates  
    % first get effective strain rates, e :
    %[etaInt,xint,yint,exx,eyy,exy,Eint,e,txx,tyy,txy]=calcStrainRatesEtaInt(CtrlVar,MUA,u,v,AGlen,n);
    [~,~,~,~,~,~,~,e,~,~,~]=calcStrainRatesEtaInt(CtrlVar,MUA,u,v,AGlen,n);
    % all these variables are are element variables defined on integration points
    % therfore if plotting on nodes, must first project these onto nodes
    eNod=ProjectFintOntoNodes(MUA,e);
    figure
    %[FigHandle,ColorbarHandel,tri]=PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,eNod,CtrlVar)    ;
    PlotNodalBasedQuantities(MUA.connectivity,MUA.coordinates,eNod,CtrlVar);
    hold on ; 
    [~,~,GLgeo]=PlotGroundingLines(CtrlVar,MUA,GF,GLgeo,xGL,yGL,'k');
    PlotMuaBoundary(CtrlVar,MUA,'b')
    title(sprintf('e at t=%-g ',CtrlVar.time)) ; xlabel('x (km)') ; ylabel('y (km)')
end


%%
if contains(plots,'-stresses-')   
    figure
    %[txzb,tyzb,txx,tyy,txy,exx,eyy,exy,e]=CalcNodalStrainRatesAndStresses(CtrlVar,MUA,AGlen,n,C,m,GF,s,b,ub,vb,ud,vd);
    [~,~,txx,tyy,txy,~,~,~,~]=CalcNodalStrainRatesAndStresses(CtrlVar,MUA,AGlen,n,C,m,GF,s,b,ub,vb,ud,vd);
    N=10;
    
    %xmin=-750e3 ; xmax=-620e3 ; ymin=1340e3 ; ymax = 1460e3 ;
    %I=find(x>xmin & x< xmax & y>ymin & y< ymax) ;
    %I=I(1:N:end);
    I=1:N:MUA.Nnodes;
    
    scale=1e-2;
    PlotTensor(x(I)/CtrlVar.PlotXYscale,y(I)/CtrlVar.PlotXYscale,txx(I),txy(I),tyy(I),scale);
    hold on
    plot(x(MUA.Boundary.Edges)/CtrlVar.PlotXYscale, y(MUA.Boundary.Edges)/CtrlVar.PlotXYscale, 'k', 'LineWidth',2) ;
    hold on ; plot(GLgeo(:,[3 4])'/CtrlVar.PlotXYscale,GLgeo(:,[5 6])'/CtrlVar.PlotXYscale,'k','LineWidth',1);
    axis equal
    axis([xmin xmax ymin ymax]/CtrlVar.PlotXYscale)
    xlabel(CtrlVar.PlotsXaxisLabel) ;
    ylabel(CtrlVar.PlotsYaxisLabel) ;    
end


%%
if contains(plots,'-sbB-')
    figure
    hold off
    AspectRatio=3; 
    ViewAndLight(1)=-40 ;  ViewAndLight(2)=20 ;
    ViewAndLight(3)=30 ;  ViewAndLight(4)=50;
    %[TRI,DT]=Plot_sbB(CtrlVar,MUA,s,b,B,TRI,DT,AspectRatio,ViewAndLight);
    Plot_sbB(CtrlVar,MUA,s,b,B,TRI,DT,AspectRatio,ViewAndLight);
end

end
