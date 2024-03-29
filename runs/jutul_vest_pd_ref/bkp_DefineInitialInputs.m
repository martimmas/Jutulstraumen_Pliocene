
function [UserVar,CtrlVar,MeshBoundaryCoordinates]=DefineInitialInputs(UserVar,CtrlVar)


%% Select the type of run by uncommenting one of the following options:
UserVar.OutputExpName='jutul_vest_lgm_TraCE-21ka';

if isempty(UserVar) || ~isfield(UserVar,'RunType')
    
%     UserVar.RunType='Inverse-MatOpt';
    % UserVar.RunType='Inverse-ConjGrad';
    % UserVar.RunType='Inverse-SteepestDesent';
    % UserVar.RunType='Inverse-ConjGrad-FixPoint';
%     UserVar.RunType='Forward-Diagnostic';
%     UserVar.RunType='Forward-Transient';
%     UserVar.RunType='TestingMeshOptions';
end

%%
% This run requires some additional input files. They are too big to be kept on Github so you
% will have to get those separately. 
%
% You can get these files on OneDrive using the link: 
%
%   https://livenorthumbriaac-my.sharepoint.com/:f:/g/personal/hilmar_gudmundsson_northumbria_ac_uk/EgrEImnkQuJNmf1GEB80VbwB1hgKNnRMscUitVpBrghjRg?e=yMZEOs
%
% Put the OneDrive folder `Interpolants' into you directory so that it can be reaced as ../Interpolants with respect to you rundirectory. 
%
%
% UserVar.IOFilesPath = 'C:/data/ua_files';
UserVar.IOFilesPath = '/proj/bolinc/users/x_marma/ua_files';
% UserVar.GeometryInterpolant=[UserVar.IOFilesPath,'/Interpolants/ant_wide/BedMachineGriddedInterpolants.mat'];
UserVar.GeometryInterpolant=[UserVar.IOFilesPath,'/Interpolants/jutul_vest/BedMachine210914_GeometryInterpolant_jutul_vest.mat'];
% UserVar.DensityInterpolant=[UserVar.IOFilesPath,'/Interpolants/ant_wide/DepthAveragedDensityGriddedInterpolant.mat'];
UserVar.SurfaceVelocityInterpolant=[UserVar.IOFilesPath,'/Interpolants/jutul_vest/Measures_SurfVelInterpolant_jutul_vest.mat'];

UserVar.CFile='FC.mat'; UserVar.AFile='FA.mat';

%%

CtrlVar.Experiment=UserVar.RunType;

switch UserVar.RunType
    
    case {'Inverse-MatOpt','Inverse-ConjGrad','Inverse-MatOpt-FixPoint','Inverse-ConjGrad-FixPoint','Inverse-SteepestDesent'}
        
        CtrlVar.InverseRun=1;
        CtrlVar.Restart=0;
        CtrlVar.Inverse.InfoLevel=1;
        CtrlVar.InfoLevelNonLinIt=0; % no information about the non-linear iterations - good to when you know all is fine
        CtrlVar.InfoLevel=0;
        UserVar.Slipperiness.ReadFromFile=0;
        UserVar.AGlen.ReadFromFile=0;
        CtrlVar.ReadInitialMesh=1;
        CtrlVar.AdaptMesh=0;
        
        CtrlVar.Inverse.Iterations=2500;
        CtrlVar.Inverse.InvertFor='logAGlenlogC' ; % {'C','logC','AGlen','logAGlen'}
        CtrlVar.Inverse.Regularize.Field=CtrlVar.Inverse.InvertFor;
        
        if contains(UserVar.RunType,'FixPoint')
            
            % FixPoint inversion is an ad-hoc method of estimating the gradient of the cost function with respect to C.
            % It can produce quite good estimates for C using just one or two inversion iterations, but then typically stagnates.
            % The FixPoint method can often be used right at the start of an inversion to get a reasonably good C estimate,
            % after which in a restart step one can switch to gradient calculation using adjoint 
            CtrlVar.Inverse.DataMisfit.GradientCalculation='FixPoint' ;
            CtrlVar.Inverse.InvertFor='logC' ;
            CtrlVar.Inverse.Iterations=1;
            CtrlVar.Inverse.Regularize.Field=CtrlVar.Inverse.InvertFor;
          
        end
        
        
    case 'Forward-Transient'
        
        CtrlVar.InverseRun=0;
        CtrlVar.TimeDependentRun=1;
        CtrlVar.Restart=0;
        CtrlVar.InfoLevelNonLinIt=1;
        UserVar.Slipperiness.ReadFromFile=1;
        UserVar.AGlen.ReadFromFile=1;
        CtrlVar.ReadInitialMesh=1;
        CtrlVar.AdaptMesh=0;
        
    case 'Forward-Diagnostic'
               
        CtrlVar.InverseRun=0;
        CtrlVar.TimeDependentRun=0;
        CtrlVar.Restart=0;
        CtrlVar.InfoLevelNonLinIt=1;
        UserVar.Slipperiness.ReadFromFile=0;
        UserVar.AGlen.ReadFromFile=0;
        CtrlVar.ReadInitialMesh=1;
        CtrlVar.AdaptMesh=0;
        
    case 'TestingMeshOptions'
        
        CtrlVar.TimeDependentRun=0;  % {0|1} if true (i.e. set to 1) then the run is a forward transient one, if not
        CtrlVar.InverseRun=0;
        CtrlVar.Restart=0;
        CtrlVar.ReadInitialMesh=1;
        CtrlVar.AdaptMesh=1;
        UserVar.Slipperiness.ReadFromFile=0;
        UserVar.AGlen.ReadFromFile=0;
        CtrlVar.AdaptMesh=1;
        CtrlVar.AdaptMeshInitial=1  ;       % remesh in first iteration (Itime=1)  even if mod(Itime,CtrlVar.AdaptMeshRunStepInterval)~=0.
        CtrlVar.AdaptMeshAndThenStop=1;    % if true, then mesh will be adapted but no further calculations performed
        % useful, for example, when trying out different remeshing options (then use CtrlVar.doAdaptMeshPlots=1 to get plots)
        CtrlVar.InfoLevelAdaptiveMeshing=10;
end


CtrlVar.dt=0.01;
CtrlVar.time=0;
CtrlVar.TotalNumberOfForwardRunSteps=1e1; 
CtrlVar.TotalTime=10;
CtrlVar.DefineOutputsDt = 2;       % every years

% Element type
CtrlVar.TriNodes=3 ;


%%
CtrlVar.doplots=0;
CtrlVar.PlotMesh=0;  
CtrlVar.PlotBCs=0;
CtrlVar.PlotXYscale=1000;
% CtrlVar.doAdaptMeshPlots=5; 
%%

UserVar.MeshOutlinePath = [UserVar.IOFilesPath,'/MeshFiles/mesh_lines_jutul_vest/footprint_jutul_vest.dat'];
CtrlVar.ReadInitialMeshFileName=[UserVar.IOFilesPath,'/MeshFiles/mesh_jutul_vest_ref_initial_dbde']; 
CtrlVar.MaxNumberOfElements=70e3;




%% Meshing 


CtrlVar.MeshRefinementMethod='explicit:local:newest vertex bisection';   
%CtrlVar.MeshRefinementMethod='explicit:local:red-green';
% CtrlVar.MeshRefinementMethod='explicit:global';   

% CtrlVar.MeshGenerator='gmsh' ; % 'mesh2d';
CtrlVar.MeshGenerator='mesh2d' ; % 'mesh2d';
CtrlVar.GmshMeshingAlgorithm=8; 
CtrlVar.MeshSizeMax=40e3;
CtrlVar.MeshSize=CtrlVar.MeshSizeMax/2;
CtrlVar.MeshSizeMin=1e3;
UserVar.MeshSizeIceShelves=CtrlVar.MeshSize/2;
% MeshBoundaryCoordinates=CreateMeshBoundaryCoordinatesForPIGandTWG(CtrlVar);
% MeshBoundaryCoordinates=CreateMeshBoundaryCoordinatesForAntartica(UserVar,CtrlVar);   
MeshBoundaryCoordinates=CreateMeshBoundaryCoordinates(CtrlVar,UserVar.MeshOutlinePath);
CtrlVar.AdaptMeshInitial=0  ;       % remesh in first iteration (Itime=1)  even if mod(Itime,CtrlVar.AdaptMeshRunStepInterval)~=0.
CtrlVar.AdaptMeshAndThenStop=0;    % if true, then mesh will be adapted but no further calculations performed
                                   % useful, for example, when trying out different remeshing options (then use CtrlVar.doAdaptMeshPlots=1 to get plots)
CtrlVar.AdaptMeshMaxIterations=5;
CtrlVar.SaveAdaptMeshFileName=[UserVar.IOFilesPath,'/MeshFiles/mesh_jutul_vest_ref_adapt']; 
CtrlVar.AdaptMeshRunStepInterval=1 ; % remesh whenever mod(Itime,CtrlVar.AdaptMeshRunStepInterval)==0


CtrlVar.MeshAdapt.GLrange=[6000 2000; 2000 1000];

I=1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Name='effective strain rates';
CtrlVar.ExplicitMeshRefinementCriteria(I).Scale=0.001;
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMin=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMax=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).p=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).InfoLevel=1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Use=true;

I=I+1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Name='lower surface gradient';
CtrlVar.ExplicitMeshRefinementCriteria(I).Scale=0.01;
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMin=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMax=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).p=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).InfoLevel=1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Use=true;

I=I+1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Name='thickness gradient';
CtrlVar.ExplicitMeshRefinementCriteria(I).Scale=0.01;
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMin=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMax=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).p=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).InfoLevel=1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Use=false;

I=I+1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Name='flotation';
CtrlVar.ExplicitMeshRefinementCriteria(I).Scale=0.0001;
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMin=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).EleMax=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).p=[];
CtrlVar.ExplicitMeshRefinementCriteria(I).InfoLevel=1;
CtrlVar.ExplicitMeshRefinementCriteria(I).Use=false;



%%
                                                        
%%  Bounds on C and AGlen
%CtrlVar.AGlenmin=1e-10; CtrlVar.AGlenmax=1e-5;
%CtrlVar.Cmin=1e-6;  CtrlVar.Cmax=1e20;        
%CtrlVar.CisElementBased=0;   
%CtrlVar.AGlenisElementBased=0;   


%% Testing adjoint parameters, start:
CtrlVar.Inverse.TestAdjoint.isTrue=0; % If true then perform a brute force calculation 
                                      % of the directional derivative of the objective function.  
CtrlVar.TestAdjointFiniteDifferenceType='central-second-order' ; % {'central','forward'}
CtrlVar.Inverse.TestAdjoint.FiniteDifferenceStepSize=1e-8 ;
CtrlVar.Inverse.TestAdjoint.iRange=[100,121] ;  % range of nodes/elements over which brute force gradient is to be calculated.
                                         % if left empty, values are calulated for every node/element within the mesh. 
                                         % If set to for example [1,10,45] values are calculated for these three
                                         % nodes/elements.
% end, testing adjoint parameters. 


if contains(UserVar.RunType,'MatOpt')
    CtrlVar.Inverse.MinimisationMethod='MatlabOptimization';
else
    CtrlVar.Inverse.MinimisationMethod='UaOptimization';
    if contains(UserVar.RunType,'ConjGrad')
        CtrlVar.Inverse.GradientUpgradeMethod='ConjGrad' ; %{'SteepestDecent','ConjGrad'}
    else
        CtrlVar.Inverse.GradientUpgradeMethod='SteepestDecent' ; %{'SteepestDecent','ConjGrad'}
    end
    
end


%%

CtrlVar.Inverse.AdjointGradientPreMultiplier='I'; % {'I','M'}


                                                    
UserVar.AddDataErrors=0;


CtrlVar.Inverse.Regularize.C.gs=1;
CtrlVar.Inverse.Regularize.C.ga=1;
CtrlVar.Inverse.Regularize.logC.ga=1;
CtrlVar.Inverse.Regularize.logC.gs=1e4 ;

CtrlVar.Inverse.Regularize.AGlen.gs=1;
CtrlVar.Inverse.Regularize.AGlen.ga=1;
CtrlVar.Inverse.Regularize.logAGlen.ga=1;
CtrlVar.Inverse.Regularize.logAGlen.gs=1e4 ;


%%
CtrlVar.ThicknessConstraints=0;
CtrlVar.ResetThicknessToMinThickness=1;  % change this later on
CtrlVar.ThickMin=10;

%%
% filename=sprintf('IR-%s-%s-Nod%i-%s-%s-Cga%f-Cgs%f-Aga%f-Ags%f-%i-%i-%s',...
%     UserVar.RunType,...
%     CtrlVar.Inverse.MinimisationMethod,...
%     CtrlVar.TriNodes,...
%     CtrlVar.Inverse.AdjointGradientPreMultiplier,...
%     CtrlVar.Inverse.DataMisfit.GradientCalculation,...
%     CtrlVar.Inverse.Regularize.logC.ga,...
%     CtrlVar.Inverse.Regularize.logC.gs,...
%     CtrlVar.Inverse.Regularize.logAGlen.ga,...
%     CtrlVar.Inverse.Regularize.logAGlen.gs,...
%     CtrlVar.CisElementBased,...
%     CtrlVar.AGlenisElementBased,...
%     CtrlVar.Inverse.InvertFor);
% filename=replace(filename,'.','k');
filename=sprintf('%s/outputs/Restart-%s-%s-%07i',UserVar.IOFilesPath,UserVar.OutputExpName,UserVar.RunType,round(100*CtrlVar.TotalTime));
CtrlVar.Inverse.NameOfRestartOutputFile=filename;
CtrlVar.Inverse.NameOfRestartInputFile=CtrlVar.Inverse.NameOfRestartOutputFile;

end
