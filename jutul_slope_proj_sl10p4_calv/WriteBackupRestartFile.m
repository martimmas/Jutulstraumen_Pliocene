function WriteBackupRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo)
%% Add to DefineOutputs.m:
%
% if isfield(UserVar,'WriteBackupRestartDt') && mod(CtrlVar.time,UserVar.WriteBackupRestartDt) == 0
%        WriteBackupRestartFile(UserVar,CtrlVar,MUA,BCs,F,GF,l,RunInfo)
% end
%
% and define UserVar.WriteBackupRestartDt in DefineInitialInputs.m!

%%

    RestartFile=sprintf('%s/outputs/RestartBkp-%s-%s-%07i',UserVar.IOFilesPath,...
                        UserVar.OutputExpName,UserVar.RunType,round(100*CtrlVar.time));    

    fprintf(CtrlVar.fidlog,' \n ################## %s %s ################### \n Writing backup restart file %s  at t=%-g \n %s \n ',CtrlVar.Experiment,datestr(now),RestartFile,CtrlVar.time);
    
    CtrlVarInRestartFile=CtrlVar;
    UserVarInRestartFile=UserVar;
    
    if isfield(RunInfo.Forward,'dtRestart')
        CtrlVarInRestartFile.dt=RunInfo.Forward.dtRestart;
    end

    
    try
        save(RestartFile,'CtrlVarInRestartFile','UserVarInRestartFile','MUA','BCs','F','GF','l','RunInfo','-v7.3');
        
        fprintf(CtrlVar.fidlog,' Writing backup restart file was successful. \n');
        
    catch exception
        fprintf(CtrlVar.fidlog,' Could not save backup restart file %s \n ',RestartFile);
        fprintf(CtrlVar.fidlog,'%s \n',exception.message);
    end
    
    