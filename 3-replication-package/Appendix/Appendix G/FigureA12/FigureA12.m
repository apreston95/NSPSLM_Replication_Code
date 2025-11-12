%KurmannSims_Part2.m
%
%
%   Kurmann, A. and E. Sims (2019). "Revisions in Utilization-Adjusted TFP and Robust Identification of News Shocks." Review of Economics and Statistics (forthcoming).
%
% Written by Andre Kurmann (Drexel University) and Eric Sims (Notre Dame
% and NBER) except for a few utilities from J.P. Lesage's econometrics
% toolbox (see acknowledgment in code).
%
%
% Use of code for research purposes is permitted as long as proper reference to source is given. 

% AP: Note that the news shock is contained in the vector 'v1'
%
%------------------------------------------------------------------------------------------------------------------------------------------------

clear all
close all

rng('default')      %reset random number generator for draws so as to always get same results

%% Step 1: parameter settings
%-------------------------------------------------------------
shock_extract = 2;      % 1: Uhlig / Francis et al max-share identification
                        % 2: Barsky-Sims TFP news identification 
                   
est_method = 1;         % 0 : OLS point estimation
                        % 1 : Bayesian estimation with Minnesota prior 
                        
sample = 1;             % 1 = quarterly macro data
                        % 2 = annual macro data (to use Alexopoulos' tech indicators for Figure 6a)                        

if est_method == 1
    prior=inf;            % for prior,  1 =  standard Minnesota prior
                        %             inf = flat prior
    ndraws = 1000;      % number of draws when computing error bands
    bound=0.16;         % lower bound of confidence interval
                        % (e.g. for bound=0.16, CI is 16% - 84%)
end
        
if sample == 1
    nlags = 4;                    % # of lags for quarterly VAR
    nhorizon = 80;                % horizon (quarters) of FEVs and IRFs to be computed
    nimp = 20;                    % quarters of FEVs and IRFs to be plotted (cannot exceed nhorizon)
    if shock_extract == 1
        KLbar=80;                 % lower horizon bound (quarters) for FEV share objective for max-share identification
        KUbar=80;                 % higher horizon bound (quarters) for FEV share objective for max-share identification
    elseif shock_extract == 2
        KLbar=0;                  % lower horizon bound (quarters) for FEV share objective for Barsky-Sims identification
        KUbar=40;                 % higher horizon bound (quarters) for FEV share objective for Barsky-Sims identification
    end
elseif sample == 2
    nlags = 2;                    % # of lags for annual VAR
    nhorizon = 20;                % horizon (years) of FEVs and IRFs to be computed
    nimp = 10;                    % years of FEVs and IRFs to be plotted (cannot exceed nhorizon)
    if shock_extract == 1
        KLbar=20;                 % lower horizon bound (years) for FEV share objective for max-share identification
        KUbar=20;                 % higher horizon bound (years) for FEV share objective for max-share identification
    elseif shock_extract == 2
        KLbar=0;                  % lower horizon bound (years) for FEV share objective for Barsky-Sims identification
        KUbar=10;                 % higher horizon bound (years) for FEV share objective for Barsky-Sims identification
    end
end

                
%% Step 2: load and transform data
%-----------------------------------------------------------------------
 
%quarterly macro data
     esty1 = 1960;     % First Year for estimation (1960:1 is first available date for UMich conf index, 1968:3 is the first date for BOS conf index)
     estq1 = 1;        % First quarter for estimation  
     esty2 = 2019;     % Last Year for estimation (2007:3 is last available date for Fernald's 2007 vintage)      
     estq2 = 4;        % Last quarter for estimation 
     [y, x, vars, lev] = dataset_quarterly_bigVAR(esty1,estq1,esty2,estq2,nlags);    
      



%% Step 3: Estimate VAR, extract shocks and VDs, IRFS 
%------------------------------------------------------

%add constant
[T,nvars]=size(y);
const = ones(T,1); 
x = [x const];
ncoeffs = nvars*nlags+1;

%ncoeffs x nvars matrix of regression coefficients
b=x\y;       %OLS coefficients: ncoeffs x nvars (more efficient way of computing least squares than b=inv(x'*x)*x'*y)
res = y - x*b;          %T x nvar matrix of residuals
vmat = (1/T)*(res'*res);
stds=sqrt(diag(vmat));

%OLS point estimates
if est_method==0; 
    if (shock_extract == 1)
        [data,ysim,v1]=uhlig(y,x,b,res,vmat,nvars,nlags,KLbar,KUbar,nhorizon);
    elseif (shock_extract == 2)
        [data,ysim,v1]=barskysims(y,x,b,res,vmat,nvars,nlags,KLbar,KUbar,nhorizon);
    end

    %data format:  rows = nimp
    %              cols = 2*nvar (FEVs and IRFs of each variable)
    %              sheets = 1
    %              4th = draws (if Bayesian) 
   
%Minnesota prior point estimates and draws
elseif est_method==1;
        
    [hm,bm] = minneprc(y,x,nlags,1,1,lev,prior); %Minnesota prior
    xxx = inv(kron(eye(nvars),x'*x) + diagrv(eye(rows(hm)),hm));
    bb = xxx * (vec(x'*y) + bm);    %stacked posterior means of regression coeffs: vec(bpost) = inv[inv(H) + kron(I,X'X)] * [vec[X'Y] + inv(H)*bprior]
    b=reshape(bb,ncoeffs,nvars);    %posterior coefficient matrix 
    res = y - x*b;                  %T x nvar matrix of residuals
    res = res - ones(T,1)*mean(res);    %%%%%%%%%
    vmat = (1/T)*(res'*res);
    
    if ndraws==0;
        if (shock_extract == 1)
            [data,ysim,v1]=uhlig(y,x,b,res,vmat,nvars,nlags,KLbar,KUbar,nhorizon);
        elseif (shock_extract == 2)    
            [data,ysim,v1]=barskysims(y,x,b,res,vmat,nvars,nlags,KLbar,KUbar,nhorizon);
        end

    else
        sxx=chol(xxx)';
        sinv=chol(inv(vmat));
        datadraws=[];
        vdraws1=[];
        randmatrix=randn(nvars,T,ndraws);
        randvec=randn(nvars*ncoeffs,ndraws);
        naccept=0;    %counter for number of accepted draws
        for j=1:ndraws;
            [bbdraw,vmatdraw] = niw(b,sxx,sinv,T,randmatrix(:,:,j),randvec(:,j)); %drawing from posterior coefficient matrix
            bdraw=reshape(bbdraw,ncoeffs,nvars);
            res = y - x*bdraw;
            res = res - ones(T,1)*mean(res);    %%%%%%%%%
            vmatdraw = res'*res/T; 
            %compute results for each draw
            if (shock_extract == 1)
                [FEVIRFs,ysim,v1]=uhlig(y,x,bdraw,res,vmatdraw,nvars,nlags,KLbar,KUbar,nhorizon);
            elseif (shock_extract == 2)              
                [FEVIRFs,ysim,v1]=barskysims(y,x,bdraw,res,vmatdraw,nvars,nlags,KLbar,KUbar,nhorizon);
            end 
            naccept=naccept+1;
            datadraws(:,:,:,naccept)=FEVIRFs;
            ysimdraws(:,:,naccept)=ysim;
            vdraws1(:,naccept)=v1;
        end
            
        %compute median and coverage bands for each result
        [temp,temp,temp,J]=size(datadraws);
        datadraws=sort(datadraws,4);
        low=round(J*bound); high=round(J*(1-bound));
        datalow=datadraws(:,:,:,low);
        datamed=median(datadraws,4);
        datahigh=datadraws(:,:,:,high);
        data(:,:,:,1)=datamed; data(:,:,:,2)=datalow; data(:,:,:,3)=datahigh;
        
        [temp,temp,J]=size(ysimdraws);
        ysimdraws=sort(ysimdraws,3);
        ysimlow=ysimdraws(:,:,low);
        ysimmed=median(ysimdraws,3);
        ysimhigh=ysimdraws(:,:,high);
        ysim(:,:,1)=ysimmed; ysim(:,:,2)=ysimlow; ysim(:,:,3)=ysimhigh;
        
        vmed1 = median(vdraws1,2);   %median values of shock
    end
end



%% Step 4: Report results
%----------------------------------------------------------------------


% plot IRFs to extracted shock (results for Figures 2, 4, 5, 6)
plotIRF(data,vars,nimp)

print('FigureA12','-dpng');


