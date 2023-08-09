
function [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar)


%% Select the type of run by uncommenting one of the following options:
UserVar.OutputExpName='jutul_vest_pd_ref';

if isempty(UserVar) || ~isfield(UserVar,'RunType')    
   UserVar.RunType='Forward-Transient';
end

% UserVar.IOFilesPath = '/proj/bolinc/users/x_marma/ua_files';
UserVar.IOFilesPath = '/mnt/3TB/data/ua_files';

% Domain
CtrlVar.ReadInitialMesh=1;           % default 1
UserVar.MeshOutlinePath = [UserVar.IOFilesPath,'/MeshFiles/mesh_lines_jutul_vest/footprint_jutul_vest.dat'];
CtrlVar.ReadInitialMeshFileName=[UserVar.IOFilesPath,'/MeshFiles/mesh_jutul_vest_ref_initial_dbde.mat'];
CtrlVar.Parallel.uvhAssembly.spmd.isOn=0; %for parallel runs

% Boundaries
UserVar.GeometryInterpolant=[UserVar.IOFilesPath,'/Interpolants/jutul_vest/BedMachine210914_GeometryInterpolant_jutul_vest.mat'];
UserVar.SurfaceVelocityInterpolant=[UserVar.IOFilesPath,'/Interpolants/jutul_vest/Measures_SurfVelInterpolant_jutul_vest.mat'];

UserVar.CFile=[UserVar.IOFilesPath,'/InputAC/FC06_bma_jutul_vest.mat']; 
UserVar.AFile=[UserVar.IOFilesPath,'/InputAC/FA_bma_jutul_vest.mat'];
CtrlVar.SlidingLaw='Weertman'; %{'Weertman'|'Umbi'|'Cornford'} default: Weertman


% Forcings
UserVar.SMB.File=[UserVar.IOFilesPath,'/Interpolants/ant_wide/FasRACMO.mat'];

UserVar.SMB.Tclim_file=[UserVar.IOFilesPath,'/Interpolants/ant_wide/lgm/FTc_TraCE-21ka.mat'];
UserVar.SMB.Pclim_file=[UserVar.IOFilesPath,'/Interpolants/ant_wide/lgm/FPc_TraCE-21ka.mat'];
UserVar.SMB.s_ref_file=[UserVar.IOFilesPath,'/Interpolants/ant_wide/lgm/Fsr_TraCE-21ka.mat'];


if ~isfile(UserVar.GeometryInterpolant) || ~isfile(UserVar.SurfaceVelocityInterpolant)     
     fprintf('\n This run requires the additional input files: \n %s \n %s  \n \n',UserVar.GeometryInterpolant,UserVar.SurfaceVelocityInterpolant)          
end

run ../load_common_run_parameters;
CtrlVar.doplots=0;
% CtrlVar.ThickMin=10;                    % default was 50

CtrlVar.MeshSizeMax=40e3;
CtrlVar.MeshSize=CtrlVar.MeshSizeMax/2;
CtrlVar.MeshSizeMin=5e2;
UserVar.MeshSizeIceShelves=5e3;
CtrlVar.MeshAdapt.GLrange=[6000 2000; 2000 1000];


%% Duration

CtrlVar.dt=0.01;
CtrlVar.time=0;
CtrlVar.TotalTime=100;
CtrlVar.TotalNumberOfForwardRunSteps=10e12;

CtrlVar.DefineOutputsDt = 10;       % every years

%% Run-type specific
CtrlVar.Experiment=UserVar.RunType;

        
CtrlVar.InverseRun=0;
CtrlVar.TimeDependentRun = 1;
CtrlVar.Restart          = 1;
CtrlVar.ResetTime        = 0;

UserVar.SeaLevelStart     = 0; % PD sea level
UserVar.SeaLevelEnd       = 0; % LGM sea level (Lambeck et al., 2014; PNAS)
UserVar.SeaLevelStartTime = 0;   % year
UserVar.SeaLevelEndTime   = 3000; % year       

UserVar.applyDeltaBedrock = 0;
UserVar.deltaBedrockInit  = 3000;
UserVar.deltaBedrockEnd   = 4000;
UserVar.deltaBedrockFile=[UserVar.IOFilesPath,'/Interpolants/ant_wide/ICE-6G_dzl_lgm.mat'];
CtrlVar.GeometricalVarsDefinedEachTransienRunStepByDefineGeometry="SB";        

CtrlVar.ATSdtMax=1.0;

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

UserVar.SMB.TimeInit=4000;
UserVar.SMB.TimeEnd=5000;
UserVar.SMB.ComputeFromTS=1;

CtrlVar.MassBalanceGeometryFeedback=3;
UserVar.Qbm.Mode='PICO'; % '{'ReadFromFile'|'Quadratic'|'PICO'}
UserVar.Qbm.GLenhance=0;

if CtrlVar.Restart==1
%     CtrlVar.NameOfRestartFiletoRead=[UserVar.IOFilesPath,'/outputs/Restart-jutul_vest_pd_ref-Forward-Transient-0001000'];
    CtrlVar.NameOfRestartFiletoRead=[UserVar.IOFilesPath,'/outputs/Restart-jutul_vest_pd_ref-Forward-Transient-0005000'];
    CtrlVar.ReadInitialMesh=0;
    CtrlVar.AdaptMesh=1;
    CtrlVar.AdaptMeshInitial=0;
    %CtrlVar.AdaptMeshTimeInterval=20;            
end

fname_restart=sprintf('%s/outputs/Restart-%s-%s-%07i',UserVar.IOFilesPath,UserVar.OutputExpName,UserVar.RunType,round(100*CtrlVar.TotalTime));
CtrlVar.NameOfRestartFiletoWrite=fname_restart;

MeshBoundaryCoordinates=CreateMeshBoundaryCoordinates(CtrlVar,UserVar.MeshOutlinePath);
%%
UserVar.AddDataErrors=0;

%% Convergence tolerances

% % The non-linear uvh/uv loops are considered to have converged if:
% %
% %  1) Work and Force tolerances are both less than: 
% CtrlVar.uvhDesiredWorkAndForceTolerances=[1e4 1e-5];
% % and, furthermore, at least one of Work and Force tolerances are less than:
% CtrlVar.uvhDesiredWorkOrForceTolerances=[1 1e-10];
% 
% % 2) If the step length in the backtracking becomes smaller than
% CtrlVar.uvhExitBackTrackingStepLength=1e-2;
% % while at the same time these Work and Force tolerances also fullfilled:
% CtrlVar.uvhAcceptableWorkAndForceTolerances=[inf 1e-3];
% CtrlVar.uvhAcceptableWorkOrForceTolerances=[1 1e-4];
% 
% 
% CtrlVar.uvDesiredWorkAndForceTolerances=[1e4 1e-5];
% CtrlVar.uvDesiredWorkOrForceTolerances=[1 1e-10];
% CtrlVar.uvExitBackTrackingStepLength=1e-2;
% CtrlVar.uvAcceptableWorkAndForceTolerances=[inf 1e-3];
% CtrlVar.uvAcceptableWorkOrForceTolerances=[1 1e-4];
% 
% CtrlVar.hDesiredWorkAndForceTolerances=[1e4 1e-5];
% CtrlVar.hDesiredWorkOrForceTolerances=[1 1e-10];
% CtrlVar.hExitBackTrackingStepLength=1e-2;
% CtrlVar.hAcceptableWorkAndForceTolerances=[inf 1e-3];
% CtrlVar.hAcceptableWorkOrForceTolerances=[1 1e-4];
% CtrlVar.hSolverMaxIterations=50;


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
