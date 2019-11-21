%% doASD_advanced5_asFunction
% Function: performs an active subspace decomposition and delivers active
% variables and acrivity scores as outputs.
% Main program of the active subspace decomposition V5
%
% Created by: Daniel Erdal, Uni Tübingen, Nov. 2019
% Last edited by: N/A
% Input:
% R         --> matrix nrX x nrParam with all paramater sets
% obs       --> vector nrX x 1 with the corresponding observations
% asdSetup  --> structure with acive subspace settings (see below)
% plotSetup --> structure with plotting settings (see below)
%
% Output:   
% aSCtot    --> activity score
% ASS       --> active subspace decomposition and surfaces
% stats     --> resulting statistics (NRMSE, validation, etc)
% asdSetup  --> as in input
% plotSetup --> as in input
%
% Sources of information:
% Constantine, P. G., Dow, E., & Wang, Q.(2014).Active Subspace Methods in Theory
%   and Practice: Applications to Kriging Surfaces.SIAM J. Sci. Comput.,
%   36(4), A1500–A1524.
% Constantine, P. G., & Diaz, P.  (2017).  Global sensitivity metrics from active subspaces.
%   Reliab. Eng. Syst. Saf.,162 (January), 1–13. doi: 10.1016/j.ress.2017.01.013
% Erdal, D., & Cirpka, O. A.    (2019). Global sensitivity analysis and adaptive stochastic
%   sampling of a subsurface-flow model using active subspaces. Hydrol. Earth Syst. Sci.,
%   23(9), 3787–3805. doi: 10.5194/hess-23-3787-2019
%
% alterations from version 4: better comments and removal of unused stuff.
% Obs: also input structure is here completely different


function [activityScore, ASS, stats, asdSetup,plotSetup]=doASD_advanced5_asFunction(R,obs,asdSetup,plotSetup) 


% asdSetup
% .gradientModel --  1 - 1st order, 2 - \todo, 3 - 3rd order polynomial,
%       4 - numerical, 5 - precalculated / 1
% .doSubSampleR -- use only the nearest neighbours to compute the gradient
%       model / false
% .nrSamplesToUse -- how many samples to use for the subsampling. Either
%       input as interger (absolute value) or ratio (beteween 0-1) / NaN
% .numericalGrad -- structure for numerical gradient / NaN
%       .model -- function handle for the model evaluation / NaN
%       .data -- input for the model
%       .eps   -- numerical epsilon / 1e-10
%       .fbc -- 1 -forward, 2- backward, 3- central difference / 1
% .Cprecalc -- precalculated C / NaN
% .valData -- validation data for validating the / NaN
%       .R -- R for validation
%       .obs -- corresponding observation
% .skipCheck -- if true, skip the pre-check (for e.g. samplers)


% plotSetup
% .doPlot -- show plot or not / false
% .figNum -- figure number to use / 100
% .plotStyle3D -- 1 - with plane, 2 - scatter3D / 2
% .xlbl -- lable for the parameters
% .title -- title to use for the figure
% .noPrinting -- if true nothing is printed throughout the run


if nargin == 1; error('This program needs at least two inputs'); end
if nargin < 4; plotSetup=[]; end
if nargin < 3; asdSetup=[]; end

% check the input and assign defaults where needed
if ~(isfield(asdSetup,'skipCheck') && asdSetup.skipCheck)
    [R,obs,asdSetup,plotSetup]=initializeASDv5(R,obs,asdSetup,plotSetup);
end
nrParam=size(R,2);


%% Get the gradients

C = computeGradients(R,obs,asdSetup);

%% Perform the Active Subspace Decomposition

% compute the eigenvector decomposition
% and order them in decreasing eigen values
[VX,DX] = eig(C);
[~,ix2]=sort(abs(diag(DX)),'descend');
V=VX(:,ix2);
d=diag(DX);
D=diag(d(ix2));
% Get the two most influential eigenvectors and their active variables
W1=V(:,1:2);
av=W1'*R';

%% Compute the activity score
if asdSetup.gradientModel==1 && ~asdSetup.doSubSampleR
    % use only one dimension
    jT=1;
else
    % use max 10 dimensions
    jT=min(10,size(D,1));
end

aSc=zeros(nrParam,jT);
for j=1:jT
    for i=1:nrParam
        aSc(i,j)=D(j,j)*(V(i,j)^2);
    end
end
aSc=cumsum(aSc,2);
aSCtot(:,1)=aSc(:,end);


%% Now for some plotting

if plotSetup.doPlot
    figure(plotSetup.figNum)
    set(plotSetup.figNum,'units','normalized','outerposition',[0.03,0,0.98,0.99])
    clf
    
    %% plot the actual eignevectors
    subplot(2,2,1);
    plot(1:nrParam,W1,'-o',[0.8,nrParam+0.2],[0,0],'--k',1:nrParam,abs(W1),'--x')
    axis tight; hold on;
    set(gca,'XTick',1:nrParam)
    set(gca,'XTickLabel',plotSetup.xlbl)
    set(gca,'XTickLabelRotation',45)
    set(gca,'FontSize',14)
    grid on
    ylabel('Weight')
    title(plotSetup.title)
    
    %% plot the activity score of the last
    subplot(2,2,3);
    hold on
    
    [~,ixASC]=sort(-aSc(:,end));
    plot(aSc(ixASC,1:end-1),'-')
    plot(aSc(ixASC,end),'-o')
    legend
    set(gca,'XTick',1:nrParam)
    set(gca,'XTickLabel',plotSetup.xlbl(ixASC))
    set(gca,'XTickLabelRotation',45)
    set(gca,'FontSize',14)
    grid on
    axis tight;
    set(gca,'XLim',[1,7])
    ylabel('Activity score (-)')
    tmp=d(ix2);
    tmp1=['1st: ' num2str(round(tmp(1),3,'significant')) ' --> ' num2str(100*round(tmp(2:5)'./tmp(1:4)',3,'significant')) ' %'];
    tmp2=round(-diff(tmp(1:5)),3,'significant')';
    title({tmp1;['Diff: ' num2str(tmp2)]})
    
    %% plot 1-D version of the active subpsace
    subplot(2,2,2);
    plot(av(1,:),obs,'o')
    
    % Fit: 'untitled fit 2'.
    [xData, yData] = prepareCurveData( av(1,:)', obs );
    % Set up fittype and options.
    ft = fittype( 'poly1' );
    % Fit model to data.
    [fitresult{1}, gof(1)] = fit( xData, yData, ft );
    % Set up fittype and options.
    ft = fittype( 'poly2' );
    % Fit model to data.
    [fitresult{2}, gof(2)] = fit( xData, yData, ft );
    ft = fittype( 'poly3' );
    % Fit model to data.
    [fitresult{3}, gof(3)] = fit( xData, yData, ft );
    % Plot fit with data.
    plot(xData, yData,'o');
    axis tight; a=axis;
    hold on
    h1=plot( fitresult{1},'r');
    h2=plot( fitresult{2},'g');
    h3=plot( fitresult{3},'m');
    set(h1,'LineWidth',2)
    set(h2,'LineWidth',2)
    set(h3,'LineWidth',2)
    ylabel('Observation')
    xlabel('Active variable')
    
    legend('obs vs. av',['Linear: R^2=' num2str(gof(1).rsquare(1))],['Qubic: R^2=' num2str(gof(2).rsquare(1))],['Quadatic: R^2=' num2str(gof(3).rsquare(1))], 'Location', 'NorthEast' );
    grid on
    axis tight
    set(gca,'FontSize',14)
    title(['Sample size: ' num2str(size(R,1))])
    stats.RsquareMax1=gof(3).rsquare(1);
    
    %% Plot the 3-D version of the active subspace if it exists
    fitresult3D=[];
    if jT~=1
        
        subplot(2,2,4);
        
        plot3(av(1,:),av(2,:),obs,'o')
        av1=av(1,:);
        av2=av(2,:);
        hold on
        a=axis;
        % Fit: 'untitled fit 1'.
        [xData, yData, zData] = prepareSurfaceData( av1', av2', obs );
        
        % Set up fittype and options.
        ft = fittype( 'poly33' );
        % Fit model to data.
        
        if plotSetup.plotStyle3D==1
            plot3(av(1,:),av(2,:),obs,'o')
            
            hold on
            a=axis;
            
            % Fit model to data.
            fitresult3D = fit( [xData, yData], zData, ft );
            
            % Plot fit with data.
            h = plot( fitresult3D);
            set(h,'FaceAlpha',0.1)
            
            % Label axes
            axis(a)
            view( -10.3, 38.6 );
            
        elseif plotSetup.plotStyle3D==2
            %
            scatter3(av(1,:),av(2,:),obs,50,obs,'filled')
            view(2)
            colorbar
            axis tight
            
        end
        
        % Label axes
        xlabel( 'av1', 'Interpreter', 'none' );
        ylabel( 'av2', 'Interpreter', 'none' );
        zlabel( 'obs', 'Interpreter', 'none' );
        grid on
        axis(a)
        view( -10.3, 38.6 );
        
        grid on
        set(gca,'FontSize',14)
        
        r=zeros(3,3);
        for i=1:3
            for j=1:3
                ft = fittype( ['poly' num2str(i) num2str(j)] );
                [fitresult3D, gof3D] = fit( [xData, yData], zData, ft );
                r(i,j)=gof3D.rsquare;
            end
        end
        stats.RsquareMax2=r(3,3);
        
        title(num2str(r))
        
       
        pred=fitresult3D(av(1,:), av(2,:));
        stats.NRMSE=sqrt(mean((((pred'-obs)).^2)./(std(obs).^2)));
        
    end
else
    % just get the highest R^2 parameter
    
    [xData, yData] = prepareCurveData( av(1,:)', obs );
    % Set up fittype and options.
    ft = fittype( 'poly3' );
    % Fit model to data.
    [fitresult{3}, gof] = fit( xData, yData, ft );
    stats.RsquareMax1=gof.rsquare;
    fitresult3D=[];
    if ~(asdSetup.gradientModel == 1 && ~asdSetup.doSubSampleR)
        [xData, yData, zData] = prepareSurfaceData( av(1,:), av(2,:), obs' );
        ft = fittype( 'poly33' );
        [fitresult3D, gof3D] = fit( [xData, yData], zData, ft );
        stats.RsquareMax2=gof3D.rsquare;        
    end
end

% make predictions and check how good they are
pred=fitresult{3}(av(1,:)');
stats.NRMSE1=sqrt(mean((((pred-obs)).^2)./(std(obs).^2)));
if ~(asdSetup.gradientModel == 1 && ~asdSetup.doSubSampleR)
    pred=fitresult3D(av(1,:), av(2,:));
    stats.NRMSE2=sqrt(mean((((pred'-obs)).^2)./(std(obs).^2)));
end


if isfield(asdSetup,'valData') && ~isempty(asdSetup.valData)
    stats.resVal1=sqrt(mean((fitresult{3}((W1'*asdSetup.valData.R')')-asdSetup.valData.obs).^2));
    if ~(asdSetup.gradientModel == 1 && ~asdSetup.doSubSampleR)
        stats.resVal2=sqrt(mean((fitresult3D((W1'*asdSetup.valData.R')')-asdSetup.valData.obs).^2));
    end
end


%% Sort out the output

activityScore=aSCtot;
ASS.fitresult3D=fitresult3D;
ASS.fitresult=fitresult;
ASS.W1=W1;



