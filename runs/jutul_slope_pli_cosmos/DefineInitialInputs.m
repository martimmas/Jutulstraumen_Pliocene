
function [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar)


%% Select the type of run by uncommenting one of the following options:
UserVar.OutputExpName='jutul_slope_pli_cosmos_fc10_pax';

if isempty(UserVar) || ~isfield(UserVar,'RunType')    
   UserVar.RunType='Forward-Transient';
end

% UserVar.IOFilesPath = '/proj/bolinc/users/x_marma/ua_files';
% UserVar.IOFilesPath = 'C:/Users/mama9638/Documents/MATLAB/ua_files';
UserVar.IOFilesPath = '/mnt/3TB/data/ua_files';

% Domain
CtrlVar.ReadInitialMesh=1;           % default 1
UserVar.MeshOutlinePath = [UserVar.IOFilesPath,'/MeshFiles/mesh_lines/footprint_jutul_slope.dat'];
CtrlVar.ReadInitialMeshFileName=[UserVar.IOFilesPath,'/MeshFiles/mesh_jutul_slope_ref_initial_dbde.mat'];
CtrlVar.Parallel.uvhAssembly.spmd.isOn=1; %for parallel runs

% Boundaries
UserVar.GeometryInterpolant=[UserVar.IOFilesPath,'/Interpolants/Jutulstraumen/slope_extent/BedMachine210914_GeometryInterpolant_jutul_slope.mat'];
UserVar.BedrockInterpolant=[UserVar.IOFilesPath,'/Interpolants/ant_wide/palaeo_topo/Paxman_plio_med_3pt5ma_5km.mat'];
UserVar.SurfaceVelocityInterpolant=[UserVar.IOFilesPath,'/Interpolants/Jutulstraumen/slope_extent/Measures_SurfVelInterpolant_jutul_slope.mat'];

UserVar.CFile=[UserVar.IOFilesPath,'/InputAC/FC10_bma_jutul_dhdte10.mat']; 
UserVar.AFile=[UserVar.IOFilesPath,'/InputAC/FA_bma_jutul_dhdte10.mat'];
CtrlVar.SlidingLaw='Weertman'; %{'Weertman'|'Umbi'|'Cornford'} default: Weertman


% Forcings
UserVar.SMB.File=[UserVar.IOFilesPath,'/Interpolants/Jutulstraumen/slope_extent/racmo_smb_jutul_slope.mat'];

UserVar.SMB.Tclim_file=[UserVar.IOFilesPath,'/Interpolants/ant_wide/pliocene/FTc_COSMOS.mat'];
UserVar.SMB.Pclim_file=[UserVar.IOFilesPath,'/Interpolants/ant_wide/pliocene/FPc_COSMOS.mat'];
UserVar.SMB.s_ref_file=[UserVar.IOFilesPath,'/Interpolants/ant_wide/pliocene/Fsr_COSMOS.mat'];


if ~isfile(UserVar.GeometryInterpolant) || ~isfile(UserVar.SurfaceVelocityInterpolant)     
     fprintf('\n This run requires the additional input files: \n %s \n %s  \n \n',UserVar.GeometryInterpolant,UserVar.SurfaceVelocityInterpolant)          
end

run ../load_common_run_parameters;
CtrlVar.doplots=0;
CtrlVar.ThickMin=10;
CtrlVar.MeshSizeMin=2000;
CtrlVar.MeshAdapt.GLrange=[6000 2000];

%% Duration

CtrlVar.dt=0.1;
CtrlVar.time=0;
CtrlVar.TotalTime=25000;
CtrlVar.TotalNumberOfForwardRunSteps=10e12;

CtrlVar.DefineOutputsDt = 100;       % every years
UserVar.WriteBackupRestartDt = 1000;
%% Run-type specific
CtrlVar.Experiment=UserVar.RunType;

        
CtrlVar.InverseRun=0;
CtrlVar.TimeDependentRun = 1;
CtrlVar.Restart          = 1;
CtrlVar.ResetTime        = 0;

UserVar.SeaLevelStart     = 25; % PD sea level
% UserVar.SeaLevelEnd       = 25; % Pliocene sea level (deBoer et al., 2014; TC after Dolan et al., 2012; GMD) - PLISMIP-ANT setup
% UserVar.SeaLevelStartTime = 1000;   % year
% UserVar.SeaLevelEndTime   = 4000; % year       

UserVar.applyDeltaBedrock = 0;
% UserVar.deltaBedrockInit  = 4000;
% UserVar.deltaBedrockEnd   = 5000;
% UserVar.deltaBedrockFile=[UserVar.IOFilesPath,'/Interpolants/ant_wide/pliocene/plismip-ant_dzl.mat'];
CtrlVar.GeometricalVarsDefinedEachTransienRunStepByDefineGeometry="SB";        

CtrlVar.ATSdtMax=0.5;
CtrlVar.ATSdtMin=0.05;

CtrlVar.InfoLevelNonLinIt=1;
UserVar.Slipperiness.ReadFromFile=1; % default 1 (if no CFile provided, reads C-Estimate.mat)
UserVar.AGlen.ReadFromFile=1;        % default 0? (if no AFile provided, reads AGlen-Estimate.mat)
UserVar.Temp.ReadFromFile=0;

CtrlVar.ReadInitialMesh=1;           
CtrlVar.AdaptMeshInitial=1;
CtrlVar.AdaptMesh=1;                 
CtrlVar.AdaptMeshRunStepInterval=0; CtrlVar.AdaptMeshTimeInterval=20;
CtrlVar.AdaptMeshMaxIterations=5;

UserVar.SMB.ReadFromFile=1;

UserVar.SMB.TimeInit=0;
UserVar.SMB.TimeEnd=1000;
UserVar.SMB.ComputeFromTS=1;

UserVar.Qbm.Calving=1;        
UserVar.CalvingThreshold = 200;
CtrlVar.MassBalanceGeometryFeedback=3;
UserVar.Qbm.Mode  = 'Quadratic'; % '{'ReadFromFile'|'Quadratic'|'PICO'}
UserVar.Qbm.Temp0 = 1.2; % avg. sub-ice shelf response in de Boer et al. (2015; TC) - PLISMIP-ANT

if CtrlVar.Restart==1
%     CtrlVar.NameOfRestartFiletoRead=[UserVar.IOFilesPath,'/outputs/Restart-jutul_slope_pli_cosmos_fc10_pax-Forward-Transient-1500000'];
%     CtrlVar.NameOfRestartFiletoRead=[UserVar.IOFilesPath,'/outputs/Restart-jutul_slope_pli_cosmos_fc10_pax-Forward-Transient-2000000'];
    CtrlVar.NameOfRestartFiletoRead=[UserVar.IOFilesPath,'/outputs/Restart-jutul_slope_pli_cosmos_fc10_pax-Forward-Transient-2500000'];
    CtrlVar.ReadInitialMesh=0;
    CtrlVar.AdaptMesh=1;
    CtrlVar.AdaptMeshInitial=0;
    CtrlVar.AdaptMeshTimeInterval=500;            
end

fname_restart=sprintf('%s/outputs/Restart-%s-%s-%07i',UserVar.IOFilesPath,UserVar.OutputExpName,UserVar.RunType,round(100*CtrlVar.TotalTime));
CtrlVar.NameOfRestartFiletoWrite=fname_restart;

MeshBoundaryCoordinates=CreateMeshBoundaryCoordinates(CtrlVar,UserVar.MeshOutlinePath);
%%
UserVar.AddDataErrors=0;

%% Convergence tolerances

% The non-linear uvh/uv loops are considered to have converged if:
%
%  1) Work and Force tolerances are both less than: 
CtrlVar.uvhDesiredWorkAndForceTolerances=[1e4 1e-5];
% and, furthermore, at least one of Work and Force tolerances are less than:
CtrlVar.uvhDesiredWorkOrForceTolerances=[1 1e-10];

% 2) If the step length in the backtracking becomes smaller than
CtrlVar.uvhExitBackTrackingStepLength=1e-2;
% while at the same time these Work and Force tolerances also fullfilled:
CtrlVar.uvhAcceptableWorkAndForceTolerances=[inf 1e-3];
CtrlVar.uvhAcceptableWorkOrForceTolerances=[1 1e-4];


CtrlVar.uvDesiredWorkAndForceTolerances=[1e4 1e-5];
CtrlVar.uvDesiredWorkOrForceTolerances=[1 1e-10];
CtrlVar.uvExitBackTrackingStepLength=1e-2;
CtrlVar.uvAcceptableWorkAndForceTolerances=[inf 1e-3];
CtrlVar.uvAcceptableWorkOrForceTolerances=[1 1e-4];

CtrlVar.hDesiredWorkAndForceTolerances=[1e4 1e-5];
CtrlVar.hDesiredWorkOrForceTolerances=[1 1e-10];
CtrlVar.hExitBackTrackingStepLength=1e-2;
CtrlVar.hAcceptableWorkAndForceTolerances=[inf 1e-3];
CtrlVar.hAcceptableWorkOrForceTolerances=[1 1e-4];
CtrlVar.hSolverMaxIterations=50;

%%
filename=sprintf('IR-%s-%s-Nod%i-%s-%s-Cga%f-Cgs%f-Aga%f-Ags%f-%i-%i-%s',...
    UserVar.RunType,...
    CtrlVar.Inverse.MinimisationMethod,...
    CtrlVar.TriNodes,...
    CtrlVar.Inverse.AdjointGradientPreMultiplier,...
    CtrlVar.Inverse.DataMisfit.GradientCalculation,...
    CtrlVar.Inverse.Regularize.logC.ga,...
    CtrlVar.Inverse.Regularize.logC.gs,...
    CtrlVar.Inverse.Regularize.logAGlen.ga,...
    CtrlVar.Inverse.Regularize.logAGlen.gs,...
    CtrlVar.CisElementBased,...
    CtrlVar.AGlenisElementBased,...
    CtrlVar.Inverse.InvertFor);
filename=replace(filename,'.','k');
CtrlVar.Inverse.NameOfRestartOutputFile=filename;
CtrlVar.Inverse.NameOfRestartInputFile=CtrlVar.Inverse.NameOfRestartOutputFile;

end
