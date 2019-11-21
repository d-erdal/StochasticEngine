
% Function: Read the water balance file from a HydroGeoSphere job on a unix-cluster
% Part of: Stochastic sampling scheme version 4
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A

% Input: 
% modelSetup --> ModelSetup stucture containing all model information
% runNR      --> unique running number for each simulation (c in MASTER)

% Output:
% obs        --> vector containing the observations
% crashed    --> flag for a crash if true


% PLEASE NOTE: this function is specfic for setting up the Käsbach Lite
% model with the 32 unknown parameters used for testing the Stochastic
% Engine. Therefore the comments below are sparse and for any other model
% the full of this function would need to be rewritten to fit the
% inidividual needs!



% [obs,crashed]=modelSetup.getObs(modelSetup,currentRuns(i,2));
function [obs,crashed]=getObsHGS_KBlite_V1(modelSetup,runNR)
crashed=false;
SSflow=nan(1,15);
obs=0;

filename=[modelSetup.data.runDir num2str(runNR) '/KaesBachLite_v1o.water_balance.dat'];
fileID = fopen(filename,'r');
a=textscan(fileID,'%s');
fclose(fileID);
%% THIS SHOULD BE MANUALLY UPDATED
nrC=22;% 11 for StT1 and 12 for StT2-StT4, 22 for StT5
mtrx=zeros(ceil(size(a{1},1)/nrC),nrC)';

runNR=1;

strt=16;
for ii=strt:size(a{1},1)    
    mtrx(runNR)=str2double(a{1}{ii});
    runNR=runNR+1;
end
mtrx=mtrx';
mtrx(sum(mtrx,2)==0,:)=[];

%%
if abs(mtrx(end,1)-1.000000000000000e+10)>1e-3
    warning(['CHECK SIMULATION FOR: ' filename])
    crashed=true;
else
    
    %"Time","topFace_LU1","topFace_LU2","topFace_LU3","bottom_head_fmf","XitTop_drain","riverDrainStretch_1","riverDrainStretch_2","riverDrainStretch_3","riverDrainStretch_4","riverDrainStretch_5","riverDrainStretch_6","riverDrainStretch_7","riverDrainStretch_8","riverDrainStretch_9","ammerSouth_HydEqHead_FT","NET1 Sources/Sinks","PM","NET2 Accumulation","ERROR (NET1-NET2)","Error rel","Error
    for iii=2:16
        % top flux
        SSflow(1,iii-1)=mtrx(end,iii);
    end
    
    % construct the unique observation vector
    rc=SSflow(9)+SSflow(10)+SSflow(11)+SSflow(8)+SSflow(7)+SSflow(12);
    rtrt=-sum(SSflow(6:14))/sum(SSflow(1:3));
    obs=[SSflow(5),rc,SSflow(11),SSflow(12),rtrt,rtrt];
    
    
end












