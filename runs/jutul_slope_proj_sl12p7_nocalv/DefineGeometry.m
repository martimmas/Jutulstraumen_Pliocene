function [UserVar,s,b,S,B,alpha]=DefineGeometry(UserVar,CtrlVar,MUA,time,FieldsToBeDefined,F)

persistent FB FB2 Fs Fh FdB

x=MUA.coordinates(:,1); y=MUA.coordinates(:,2);
alpha=0 ;

if nargin<5
    FieldsToBeDefined='sbSB';
end

fprintf('DefineGeometry %s \n',FieldsToBeDefined)

if isempty(FB)
    
    %%
    %locdir=pwd;
    %AntarcticGlobalDataSets=getenv('AntarcticGlobalDataSets');
    %cd(AntarcticGlobalDataSets)
    fprintf('DefineGeometry: loading file: %-s ',UserVar.GeometryInterpolant)
    load(UserVar.GeometryInterpolant,'FB','Fh','Fs')   
    if isfield(UserVar,'ThicknessInterpolant')
        load(UserVar.ThicknessInterpolant,'Fh')        
    end    
end
if isfield(UserVar,'BedrockInterpolant') && isempty(FB2)
    FB2 = load(UserVar.BedrockInterpolant).FB;
end
fprintf(' done \n')
%cd(locdir)


if contains(FieldsToBeDefined,'S')
    if isfield(UserVar,'SeaLevelEnd') && (time >= UserVar.SeaLevelStartTime)
        sl_rate = (UserVar.SeaLevelEnd-UserVar.SeaLevelStart)/(UserVar.SeaLevelEndTime-UserVar.SeaLevelStartTime);
        z_sl = UserVar.SeaLevelStart + sl_rate*(time-UserVar.SeaLevelStartTime);
        S=ones(size(x)).*round(z_sl,2);
        
        if time >= UserVar.SeaLevelEndTime
            S=ones(size(x)).*UserVar.SeaLevelEnd;
        end
        
    elseif isfield(UserVar,'SeaLevelStart')        
        S=ones(size(x))*UserVar.SeaLevelStart;
    else
        S=x*0;
    end
else
    S=NaN;
end


if contains(FieldsToBeDefined,'s')
    s=Fs(x,y);        
else
    s=NaN;
end

b=NaN; B=NaN;

if contains(FieldsToBeDefined,'b')  || contains(FieldsToBeDefined,'B')
    
    if isfield(UserVar,'BedrockInterpolant')
        B=FB2(x,y);        
    else
        B=FB(x,y);
        h=Fh(x,y);
    end
    
    if UserVar.applyDeltaBedrock && isfield(UserVar,'deltaBedrockFile')
        if isempty(FdB)
            FdB = load(UserVar.deltaBedrockFile).F;    
        end
        % alternative 1: apply all at once
        %dB = FdB(x,y);
        %B=double(B+dB);
        
        % alternative 2: apply it progressively, moving the ice thickness
        % with it
        if time >= UserVar.deltaBedrockInit
            dB_rate = FdB(x,y)/(UserVar.deltaBedrockEnd-UserVar.deltaBedrockInit);
            if time > UserVar.deltaBedrockEnd
                time_app = UserVar.deltaBedrockEnd;
            else
                time_app = time;
            end
            dB = dB_rate * (time_app-UserVar.deltaBedrockInit);
            if isnan(s) % s needs to be updated regardless due to bedrock adjustment
                s = F.s;
            end            
            B = double(B-dB);
            [b,h,~] = Calc_bh_From_sBS(CtrlVar,MUA,s,B,S,900+zeros(MUA.Nnodes,1),1030);
            Igrd = B==b; % it is grounded where the ice base matches the bedrock elevation
            s(Igrd) = B(Igrd)+h(Igrd);
        end 
    end
            
%     if isfield(UserVar,'ThicknessInterpolant') % this does not work for floating ice!!
%         h=Fh(x,y);
%         s = B+h;        
%     end
    

    [b,~,~] = Calc_bh_From_sBS(CtrlVar,MUA,s,B,S,900+zeros(MUA.Nnodes,1),1030);
    %b=s-h;

    I=b<B; b(I)=B(I);  % make sure that interpolation errors don't create a situation where b<B
    
    % Vostok
    I=x>1000e3 & x < 1600e3 & y>-500e3 & y<-200e3  ;
    B(I)=b(I); % shift bed upwards towards lower ice surface and ignore lake
    
    % Correct for smoothing on the grounded part
    %I=y<1900e3;
    %b(I)=B(I); % shift ice base downwards to force the grounding due to low resolution        
    
end



end



