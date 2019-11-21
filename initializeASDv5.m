
% Function: initialization routine for the active subspace program
% checks input, assigns defaults and print warnings and info
% Part of: Active subspace decomposition V5
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A

% Input:
% R         --> matrix nrX x nrParam with all paramater sets
% obs       --> vector nrX x 1 with the corresponding observations
% asdSetup  --> structure with acive subspace settings (see below)
% plotSetup --> structure with plotting settings (see below)

% Output:
% same as input, but updated 


function [R,obs,asdSetup,plotSetup]=initializeASDv5(R,obs,asdSetup,plotSetup)

% asdSetup
% .gradientModel --  1 - 1st order, 2 - \todo, 3 - 3rd order polynomial, 
%       4 - numerical, 5 - precalculated / data dependent 1 or 3
% .doSubSampleR -- use only the nearest neighbours to compute the gradient
%       model / false
% .nrSamplesToUse -- how many samples to use for the subsampling. Either
%       input as interger (absolute value) or ratio (beteween 0-1) / NaN
% .numericalGrad -- structure for numerical gradient / NaN
%       .model -- function handle for the model evaluation / NaN
%       .data -- input for the model
%       .eps   -- numerical epsilon / 1e-10
% .Cprecalc -- precalculated C / NaN
% .valData -- validation data for validating the / NaN
%       .R -- R for validation
%       .obs -- corresponding observation

% plotSetup
% .doPlot -- show plot or not / false
% .figNum -- figure number to use / 100
% .plotStyle3D -- 1 - with plane, 2 - scatter3D / 2
% .xlbl -- lable for the parameters  / 1,2,3,4....
% .title -- title to use for the figure / ''
% .noPrinting -- if true nothing is printed throughout the run

% check if the sizes agree and also check if there are NaNs in the data
if size(R,1) ~= size(obs,1)
    error('Parameter matrix R and observations obs must have same length')
end

if any(isnan(obs))
   sOrg=size(obs,1);
   R(isnan(obs),:)=[];
   obs(isnan(obs),:)=[];
   
   txt=['Nans deleted: old size ' num2str(sOrg) ', new size ' num2str(size(obs,1))];    
end


if ~isfield(plotSetup,'doPlot')
    plotSetup.doPlot=false;
end
if ~isfield(plotSetup,'noPrinting')
    plotSetup.noPrinting=true;
end
if ~isfield(plotSetup,'figNum')
    plotSetup.figNum=100;
end
if ~isfield(plotSetup,'plotStyle3D')
    plotSetup.plotStyle3D=2;
end
if ~isfield(plotSetup,'xlbl') || isempty(plotSetup.xlbl)
        plotSetup.xlbl=num2cell(1:size(R,2)); 
end
if ~isfield(plotSetup,'title')
    plotSetup.title='';
end


if ~isfield(asdSetup,'gradientModel')
    asdSetup.gradientModel = 1;
end
if ~isfield(asdSetup,'doSubSampleR')
    asdSetup.doSubSampleR = false;
end
if ~isfield(asdSetup,'nrSamplesToUse')
    asdSetup.nrSamplesToUse=round(0.5*size(R,1));
elseif asdSetup.nrSamplesToUse < 1
    % we have a ratio; translate to absolute number
    asdSetup.nrSamplesToUse=round(asdSetup.nrSamplesToUse*size(R,1));    
end

if asdSetup.gradientModel == 4
    % numerical
    if ~isfield(asdSetup,'numericalGrad')
        error('To use numerical gradients, please specify field numericalGrad in adsSetup')
    end
    if ~isfield(asdSetup.numericalGrad,'data')
        error('To use numerical gradients, please specify field data in adsSetup.numericalGrad')
    end
    if ~isfield(asdSetup.numericalGrad,'model')
        asdSetup.numericalGrad.eps=1e-10;
    end
end

if asdSetup.gradientModel == 5 && ~isfield(asdSetup,'Cprecalc')
    % precalculated C
   error('To use precalculated C, please specify field Cprecalc in adsSetup')
end
   
if ~any(asdSetup.gradientModel == 1:5)
    error('Gradient model most be an integer between 1 and 5')
end
if ~plotSetup.noPrinting
    disp('Welcome to the Active Subpsace program')
    disp('This is version 5, november 2019, D. Erdal@ uni Tübingen, DE')
    disp(['Your are using gradient model ' num2str(asdSetup.gradientModel)])
    
    disp('---------asdSetup-----------')
    disp(asdSetup)
    disp('---------plotSetup-----------')
    disp(plotSetup)
    disp('-----------------------------')
    
    if exist('txt','var')~=0
        disp(txt)
        disp('-----------------------------')
    end
end
    
    
    
    
    
    
    
end

