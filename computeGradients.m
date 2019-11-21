
% Function: computes the gradients of the observation with respect to all the 
% parameters using a user selected method.
% Part of: Active subspace decomposition V5
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A
%
% Input:
% R         --> matrix nrX x nrParam with all paramater sets
% obs       --> vector nrX x 1 with the corresponding observations
% asdSetup  --> structure with acive subspace settings (see other functions)
%
% Output
% C     --> matrix of monte carlo sample of the gradient product, needed
%       for the active subspace decomposition


function C = computeGradients(R,obs,asdSetup)
nrParam=size(R,2);

if asdSetup.gradientModel == 5
    % C has been computer somewhere else (e.g. analytically with a
    % GPE-model)
    C = asdSetup.Cprecalc;
elseif asdSetup.gradientModel == 4
    % numerically compute the gradient
    grM=zeros(size(R,1),nrParam);   
    obsPred1 = asdSetup.numericalGrad.model(R,asdSetup.numericalGrad.data);
    
    for j=1:nrParam
        Rpr=R;
        Rpr(:,j)=R(:,j)+asdSetup.numericalGrad.eps;
        obsPred2 = asdSetup.numericalGrad.model(Rpr,asdSetup.numericalGrad.data);          
        grM(:,j)=(obsPred1-obsPred2)./asdSetup.numericalGrad.eps;        
    end
    C=zeros(nrParam);
    for j=1:size(R,1)
        C=C+grM(j,:)'*grM(j,:);
    end
    C=C/size(R,1);
else
    % Compute the gradients using a simple model
    
    if asdSetup.gradientModel==1
        % linear gradient model
        X=[ones(size(R,1),1),R];
    elseif asdSetup.gradientModel==2
        % use only self square terms
        X=[ones(size(R,1),1),R,zeros(size(R,1),nrParam)];
        for i=1:nrParam
            X(:,i+(nrParam+1))=R(:,i).*R(:,i);
        end
    elseif asdSetup.gradientModel==3
        % full third order polynomial
        X=[ones(size(R,1),1),R,zeros(size(R,1),nrParam^2)];
        cc=nrParam+2;
        for i=1:nrParam
            R2=repmat(R(:,i),1,nrParam);
            R2(:,1:i-1)=0;
            X(:,cc:cc+nrParam-1)=R.*R2;
            cc=cc+nrParam;
        end
        ix=X(1,:)==0;
        X(:,ix)=[];
    end
    
    if ~asdSetup.doSubSampleR
        m=(X'*X)\(X'*obs);    % SSqE
        
        mm=m;
        v=mm(2:(nrParam+1));
        mm(1:(nrParam+1))=[];
        H=zeros(nrParam);
        if asdSetup.gradientModel==3
            for i=1:nrParam
                H(i,i:end)=mm(1:(nrParam+1)-i);
                mm(1:(nrParam+1)-i)=[];
            end
        elseif asdSetup.gradientModel==2
            H=diag(mm);
        end
        
        
        H=H+H';
        x=R';
        
        C=zeros(nrParam);
        for i=1:size(R,1)
            grF=H*x(:,i)+v;
            C=C+grF*grF';
        end
        C=C/size(R,1);
        
    else
        C=zeros(nrParam);
        for i=1:size(R,1)
            % for each observation, approxiamte the gradient locally
            
            % find the clostest neighbours in parameter space
            a=sqrt(mean((R-repmat(R(i,:),size(R,1),1)).^2,2));
            [~,ixGR]=sort(a);
            
            % recompute X and m
            Xuse=X(ixGR(1:asdSetup.nrSamplesToUse),:);
            m=(Xuse'*Xuse)\(Xuse'*obs(ixGR(1:asdSetup.nrSamplesToUse)));
            
            mm=m;
            v=mm(2:(nrParam+1));
            mm(1:(nrParam+1))=[];
            H=zeros(nrParam);
            if asdSetup.gradientModel==3
                for i22=1:nrParam
                    H(i22,i22:end)=mm(1:(nrParam+1)-i22);
                    mm(1:(nrParam+1)-i22)=[];
                end
            elseif asdSetup.gradientModel==2
                H=diag(mm);
            end
            H=H+H';
            grF=H*R(i,:)'+v;
            
            C=C+grF*grF';
            
        end
        C=C/size(R,1);
    end
end