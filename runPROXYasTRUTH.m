% Function: Sets up and submitts a HydroGeoSphere job on a unix-cluster
% Part of: Stochastic sampling scheme version 4
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A

% Input: 
% newR      --> parameters to be used
% modelData --> modelSetup.data, structure containing model specfic input
%            Here: model name and data

% Output:
% pid    --> pid from the cluster: not used here
% obsX   --> observations from the model 

function [pid,obs]=runPROXYasTRUTH(newR,modelData,~)
pid=nan;

if strcmp(modelData.name,'AS')
    
    
    obs=zeros(size(newR,1),length(modelData.fitresult3D));
    for i=1:size(obs,2)     
        av=(modelData.W1{i}'*newR')';
        obs(:,i)=modelData.fitresult3D{i}(av);
    end   
    
    
elseif strcmp(modelData.name,'GPE')
    
    obs=zeros(size(newR,1),length(modelData.GPE));
    for i=1:size(obs,2)   
        o=stk_predict(modelData.GPE{i}.M_post,newR);
        obs(:,i)=o.mean;
    end 
    
else
    error('NO SUCH MODEL NAME YET')
end


