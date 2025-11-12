%uhlig.m
%
%function that extract Uhlig impulse matrices and computes variance
%decompositions and impulse responses to max-share shock
%
%based on code from Kurmann-Otrok (AER 2013)
%
%Andre Kurmann, Drexel Unviersity, last modified January 2020

function [FEVIRFs,ysim,v]=uhlig(y,x,b,res,vmat,nvars,nlags,KLbar,KUbar,nimp)

%companion matrix of demeaned VAR (don't include constant term)
M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=b(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);

% Extract impulse vectors 
%------------------------

%Lower triangular matrix s.t. vcv of fundamental shocks is I matrix = Cholesky decomp of vmat
Atilde = chol(vmat)';

%compute Rtildes for each l
Rtilde=zeros(nvars,nvars,KUbar+1);
FEV=zeros(nvars,nvars,KUbar+2);
for l=0:KUbar;
    C_l=M^l;
    Rtilde(:,:,l+1)=C_l(1:nvars,1:nvars)*Atilde;
end

%compute S
Eii=zeros(nvars,nvars); Eii(1,1)=1;     % ** target variable must be ordered first **
S=zeros(nvars,nvars);
S_l=zeros(nvars,nvars,KUbar+1);
eY=zeros(nvars,nvars*nlags);
eY(1:nvars,1:nvars)=eye(nvars);
for l=0:KUbar;
    S=S+(KUbar+1-max(KLbar,l))*Rtilde(:,:,l+1)'*Eii*Rtilde(:,:,l+1);
end

[V,D]=eig(S);               %eigenvectors have norm=1
lambda=[diag(D) seqa(1,1,nvars)];
order=sortrows(lambda,-1);      %sort eigenvalues in descending order
Vord=[];
for i=1:nvars;
    Vord=[Vord V(:,order(i,2))];    %reorder eigenvector according to ordered
end                                 %eigenvalues            
V=Vord;
                          
%impulse vectors associated with largest eigenvalue
q1=V(:,1);
alpha=Atilde*q1;   

v = q1'*inv(Atilde)*res';   %1*T series of the fundamental shock 
v = v';

% Compute impulse responses
%----------------------------------------------------------------------
U1=[alpha; zeros(nvars*nlags-nvars,1)];
    
for k=1:nimp
    Zk1(k,:)=(M^(k-1)*U1)';
end        
impulse1(:,:)=Zk1(:,1:nvars);   %nimp x nvars impulse response matrix to 1st shock

%test whether first variable in impulse response vector is positive, otherwise revert
%ATTENTION: NEED TO CORRECTLY SPECIFY DEPENDING ON IDENTIFICATION
%when this program is used to identify term structure slope or level shock; then, rotation condition should be on impact
%but if program is used to identify TFP max-share shock, then rotation condition should be imposed at longer horizon, say 40 quarters 
if impulse1(nimp/2,1,1) < 0   %in order to avert 180 degree rotations, condition on response to be positive; otherwise mirror;
        alpha = -alpha; 
        v = -v;
        impulse1 = -impulse1;
end


% Compute fraction of VD explained at different horizons
%-----------------------------------------------------------------------
sigmaks=zeros(nvars,nvars,nimp);
sigmak=Rtilde(:,:,1)*Rtilde(:,:,1)';
sigmaks(:,:,1)=sigmak;
hh1=Rtilde(:,:,1)*q1*(Rtilde(:,:,1)*q1)';
vardec1(1,:)=(diag(hh1./sigmak))';
    
for k=1:nimp-1;
        %add square of k-step ahead forecast error to build k-ahead variance-covariance 
        sigmak=sigmak+Rtilde(:,:,k+1)*Rtilde(:,:,k+1)';
        hh1=hh1+Rtilde(:,:,k+1)*q1*(Rtilde(:,:,k+1)*q1)';
        vardec1(k+1,:)=(diag(hh1./sigmak))';
end

%simulate historical time series implied by shock
[T,nvars] = size(y);
nlags = (size(x,2)-1)/nvars;
xsim = x(1,:);        
for i=1:T
    ysim(i,:) = xsim*b + alpha'*v(i);
    xsim = [ysim(i,:) xsim(1:nvars*(nlags-1)) 1];  %update vector of regressors
end

%build matrix of FEV and IRF outputs
FEVIRFs=[vardec1 impulse1];
