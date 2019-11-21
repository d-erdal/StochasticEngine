% Function: Sets up and submitts a HydroGeoSphere job on a unix-cluster
% Part of: Stochastic sampling scheme version 4
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A

% Input: 
% newR      --> parameters to be used
% modelData --> ModelSetup.data, structure containing model specfic input
%            Here: base folder and run directory   
% runNR     --> unique running number for each simulation (c in MASTER)

% Output:
% pid    --> pid from the cluster
% obsX   --> observation if directly available: not used here


% PLEASE NOTE: this function is specfic for setting up the Käsbach Lite
% model with the 32 unknown parameters used for testing the Stochastic
% Engine. Therefore the comments below are sparse and for any other model
% the full of this function would need to be rewritten to fit the
% inidividual needs!

function [pid,obsX]=runHGS_KBlite_V1(newR,modelData,runNR)
obsX=nan;

% set up the input for the stochstic run 

% create the destination folder
copyfile(modelData.HGSbaseFolder,[modelData.runDir num2str(runNR)])
cd([modelData.runDir num2str(runNR)])

% create and edit the input parameters and source layers
%% layers

bottomP=(-5+(5--5)*newR(1));
offsetMOW=(0+(100-0)*newR(2));
offsetkukm1=(-20+(20--20)*newR(3));
km1Wthickness=(5+(50-5)*newR(4));
% version 2 --> limit thickness
% km1Wthickness=(15+(50-15)*newR(4));sa
% 
[xx,yy]=meshgrid(3494405:10:3494400+10*560,5376200+10*660-5:-10:5376200);


%% simples: kukm1 offset

load('createDataSetData.mat', 'bottomFTRv2km1orderIx7KBEXT')
load('createDataSetData.mat', 'km1_CORRECTIONmatrx')

org=bottomFTRv2km1orderIx7KBEXT;
org(org<-999)=nan;
b=km1_CORRECTIONmatrx==2;
b2=[b(:,end),b(:,1:end-1)];
org(b)=org(b2);


r=org+offsetkukm1;
filename='bottomFTRv2km1orderIx7KBEXTcorr_USE.asc';
dx=10;
% Open the text file.
fileID = fopen(filename,'w');
fprintf(fileID,['ncols        ' num2str(size(r,2)) '\n']);
fprintf(fileID,['nrows        ' num2str(size(r,1)) '\n']);
fprintf(fileID,['xllcorner    ' num2str(xx(1)-dx/2) '\n']);
fprintf(fileID,['yllcorner    ' num2str(yy(end)-dx/2) '\n']);
fprintf(fileID,['cellsize     ' num2str(dx) '\n']);
fprintf(fileID,'NODATA_value  -99999\n');
r(isnan(r))=-99999;
for i=1:size(r,1)
    a=num2str(r(i,1:end));
    fprintf(fileID, [a '\n']);%, 'char');
end
fclose(fileID);
clear bottomFTRv2km1orderIx7KBEXT km1_CORRECTIONmatrx
    
%% Weathering zone thickness

load('createDataSetData.mat', 'bottomFTRv2km1WDEMm20KBEXT')
load('createDataSetData.mat', 'km1_CORRECTIONmatrx')

org=bottomFTRv2km1WDEMm20KBEXT+20;

org(org<-999)=nan;
b=km1_CORRECTIONmatrx==2;
b2=[b(:,end),b(:,1:end-1)];
org(b)=org(b2);
filename='bottomFTRv2km1WDEMKBEXTcorr_USE.asc';

r=org-km1Wthickness;
dx=10;
% OBS: y must be going from high to low!!!
% Open the text file.
fileID = fopen(filename,'w');

fprintf(fileID,['ncols        ' num2str(size(r,2)) '\n']);
fprintf(fileID,['nrows        ' num2str(size(r,1)) '\n']);
fprintf(fileID,['xllcorner    ' num2str(xx(1)-dx/2) '\n']);
fprintf(fileID,['yllcorner    ' num2str(yy(end)-dx/2) '\n']);
fprintf(fileID,['cellsize     ' num2str(dx) '\n']);
fprintf(fileID,'NODATA_value  -99999\n');
r(isnan(r))=-99999;

for i=1:size(r,1)
    a=num2str(r(i,1:end));
    fprintf(fileID, [a '\n']);%, 'char');
end

fclose(fileID);
clear bottomFTRv2km1WDEMm20KBEXT km1_CORRECTIONmatrx


%% Bottom pressure

load('createDataSetData','modeltop','bottomNodesXYZ','yMF','xMF','hydhead','topNodesXYZ')
z=-30+2.5:5:modeltop;

filename='bottomHeadFromMFnoArti_USE.lst';
% Open the text file.
fileID = fopen(filename,'w');
fprintf(fileID, [num2str(size(bottomNodesXYZ,1)) '\n']);
hSave=zeros(size(bottomNodesXYZ,1),1);
for i=1:size(bottomNodesXYZ,1)
    
    iy=find(min(abs(bottomNodesXYZ(i,2)-yMF(:,1)))==abs(bottomNodesXYZ(i,2)-yMF(:,1)));
    ix=find(min(abs(bottomNodesXYZ(i,1)-xMF(1,:)))==abs(bottomNodesXYZ(i,1)-xMF(1,:)));
    iz=find(min(abs(bottomNodesXYZ(i,3)-z))==abs(bottomNodesXYZ(i,3)-z));   
    
    hSave(i)=mean(mean(mean(hydhead(iy,ix,iz))))+bottomP; % mean here only if the choice ix/iy/iz is more than one
    
    if hSave(i) > (topNodesXYZ(i,3)-5)
        % atisian conditions are not wanted
        hSave(i)=topNodesXYZ(i,3)-5;
    end

fprintf(fileID, [num2str(hSave(i)) '\n']);

end
fclose(fileID);
clear modeltop bottomNodesXYZ yMF xMF hydhead topNodesXYZ

%% MO-west offset


load('createDataSetData.mat', 'bottom_ku_noFaults')
load('createDataSetData.mat', 'zoneMapKBL')

org=bottom_ku_noFaults;% bottomFTRv2kuorderIx6KBEXTPLANES;
org(org<-999)=nan;

filename='bottomFTRv2ku_NOFAULTSBASE_USE.asc';

r=org;
r2=org+offsetMOW;
r(zoneMapKBL==2)=r2(zoneMapKBL==2);
dx=10;
% Open the text file.
fileID = fopen(filename,'w');

fprintf(fileID,['ncols        ' num2str(size(r,2)) '\n']);
fprintf(fileID,['nrows        ' num2str(size(r,1)) '\n']);
fprintf(fileID,['xllcorner    ' num2str(xx(1)-dx/2) '\n']);
fprintf(fileID,['yllcorner    ' num2str(yy(end)-dx/2) '\n']);
fprintf(fileID,['cellsize     ' num2str(dx) '\n']);
fprintf(fileID,'NODATA_value  -99999\n');
r(isnan(r))=-99999;

for i=1:size(r,1)
    a=num2str(r(i,1:end));
    fprintf(fileID, [a '\n']);%, 'char');
end

fclose(fileID);
clear bottom_ku_noFaults zoneMapKBL


%% parameters
% load subsurface properties template
load SSP SSP

% conductivities
i4=log(0.5e-4);i5=log(1e-5); i7=log(1e-7); i8=log(1e-8); i9=log(1e-9); 

k=exp(i4+(i7-i4)*newR(5));
SSP=replace(SSP,'&&K_km1W_XY&&',num2str(k));
SSP=replace(SSP,'&&K_km1W_Z&&',num2str(k/(1+(50-1)*newR(6)))); % now a ratio
k=exp(i5+(i8-i5)*newR(7));
SSP=replace(SSP,'&&K_ku_XY&&',num2str(k));
SSP=replace(SSP,'&&K_ku_Z&&',num2str(k/(1+(50-1)*newR(8))));
k=exp(i7+(i9-i7)*newR(9));
SSP=replace(SSP,'&&K_km1_XY&&',num2str(k));
SSP=replace(SSP,'&&K_km1_Z&&',num2str(k/(1+(50-1)*newR(10))));
SSP=replace(SSP,'&&K_mo&&',num2str(exp(i5+(i7-i5)*newR(11))));
SSP=replace(SSP,'&&K_Q&&',num2str(exp(i5+(i7-i5)*newR(12))));

clear i*

% alpha
iL=0.5 ;iH=5; 

SSP=replace(SSP,'&&A_km1W&&',num2str((iL+(iH-iL)*newR(13))));
SSP=replace(SSP,'&&A_ku&&',num2str((iL+(iH-iL)*newR(14))));
SSP=replace(SSP,'&&A_km1&&',num2str((iL+(iH-iL)*newR(15))));
SSP=replace(SSP,'&&A_mo&&',num2str((iL+(iH-iL)*newR(16))));
SSP=replace(SSP,'&&A_Q&&',num2str((1+(5-1)*newR(17))));

clear i*

% n

iL=1.5 ;iH=9; 

SSP=replace(SSP,'&&n_km1W&&',num2str((iL+(iH-iL)*newR(18))));
SSP=replace(SSP,'&&n_ku&&',num2str((iL+(iH-iL)*newR(19))));
SSP=replace(SSP,'&&n_km1&&',num2str((iL+(iH-iL)*newR(20))));
SSP=replace(SSP,'&&n_mo&&',num2str((iL+(iH-iL)*newR(21))));
SSP=replace(SSP,'&&n_Q&&',num2str((1.5+(3-1.5)*newR(22))));

clear i*

% spec. store
iL=1e-6 ;iH=1e-4; 

SSP=replace(SSP,'&&ss_km1W&&',num2str((iL+(iH-iL)*newR(23))));
SSP=replace(SSP,'&&ss_ku&&',num2str((iL+(iH-iL)*newR(24))));
SSP=replace(SSP,'&&ss_km1&&',num2str((iL+(iH-iL)*newR(25))));
SSP=replace(SSP,'&&ss_mo&&',num2str((iL+(iH-iL)*newR(26))));
SSP=replace(SSP,'&&ss_Q&&',num2str((iL+(iH-iL)*newR(27))));
%

fid=fopen('subsurfProp_v2.mprops','w');
for i=1:size(SSP,1)
    fprintf(fid,char(SSP(i,:)));
    fprintf(fid,newline);
end
fclose(fid);

%% Boundaries
dr=(0.005+(0.2-0.005)*newR(28));
aex=(335+(355-335)*newR(29));

d=1000*3600*24*365;
Rcrop=(100 + (150-100)*newR(30))/d;
Rgrass=(80 + (130-80)*newR(31))/d;
Rforest=(100 + (150-100)*newR(32))/d;


load boundaries.mat AmmerRiverBoundary
AmmerRiverBoundary=replace(AmmerRiverBoundary,'&&aex&&',num2str(aex));

fid=fopen('AmmerRiverBoundary.txt','w');
for i=1:size(AmmerRiverBoundary,1)
    fprintf(fid,char(AmmerRiverBoundary(i,:)));
    fprintf(fid,newline);
end
fclose(fid);
clear AmmerRiverBoundary

load boundaries.mat TopFlowBoundaries
TopFlowBoundaries=replace(TopFlowBoundaries,'&&Rcrop&&',num2str(Rcrop));
TopFlowBoundaries=replace(TopFlowBoundaries,'&&Rgrass&&',num2str(Rgrass));
TopFlowBoundaries=replace(TopFlowBoundaries,'&&Rforest&&',num2str(Rforest));

fid=fopen('TopFlowBoundaries.txt','w');
for i=1:size(TopFlowBoundaries,1)
    fprintf(fid,char(TopFlowBoundaries(i,:)));
    fprintf(fid,newline);
end
fclose(fid);
clear TopFlowBoundaries


load boundaries.mat KaesbachRiverBoundaries
KaesbachRiverBoundaries=replace(KaesbachRiverBoundaries,'&&dr&&',num2str(dr));
fid=fopen('KaesbachRiverBoundaries.txt','w');
for i=1:size(KaesbachRiverBoundaries,1)
    fprintf(fid,char(KaesbachRiverBoundaries(i,:)));
    fprintf(fid,newline);
end
fclose(fid);
clear KaesbachRiverBoundaries

save newR newR

%% submit the job

cm=['sed -i "s/&&rr&&/' num2str(runNR) '/g" "qbatchX"'];
system(cm);
[~,pid]=system('qsub qbatchX');

pid=str2double(pid);

system('rm *.mat');


cd ..
























