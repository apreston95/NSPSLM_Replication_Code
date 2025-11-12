% barskysims.m
%
% function that extract Barsky-Sims news shock and computes variance
% decompositions and impulse responses.
%
% based on code from Kurmann-Otrok (AER 2013)
%
% Andre Kurmann, Drexel Unviersity, last modified January 2020

function [FEVIRFs,ysim,v] = barskysims(y,x,b,res,vmat,nvars,nlags,KLbar,KUbar,nimp)

%companion matrix of demeaned VAR (don't include constant term)
M=zeros(nvars*nlags,nvars*nlags);			
M(1:nvars,:)=b(1:nvars*nlags,:)';
M(nvars+1:nvars*nlags,1:nvars*nlags-nvars)=eye(nvars*nlags-nvars);

% Extract impulse vectors 
%------------------------

%Lower triangular matrix s.t. vcv of fundamental shocks is I matrix = Cholesky decomp of vmat
Gtilde = chol(vmat)';

% compute Bi for each l
ei = zeros(nvars,1); 
ei(1) = 1;      %selection matrix; TFP needs to be ordered first in VAR

nhorizon = max(nimp,KUbar);
B = zeros(nvars,nvars,nhorizon);
Bi = zeros(1,nvars,nhorizon+1);
for l=0:nhorizon;
    C=M^l;
    B(:,:,l+1) = C(1:nvars,1:nvars);
    Bi(:,:,l+1) = ei'*C(1:nvars,1:nvars); % Bi are MA coeffs on how ith variable responds to different shocks l+1 periods back
                                       % Bi(:,:,l+1) is the coef on epsilon(t-l)
end

%compute V
V = zeros(nvars,nvars);
for l=0:KUbar;
    V = V + (KUbar+1-max(KLbar,l))*(Bi(:,:,l+1)*Gtilde)'*(Bi(:,:,l+1)*Gtilde);
end
VHat = V(2:nvars,2:nvars);

[eigenVector,eigenValue]=eig(VHat);               %eigenvectors have norm=1
lambda=[diag(eigenValue) seqa(1,1,nvars-1)];
order=sortrows(lambda,-1);     %sort eigenvalues in descending order
Vord=[];
for i=1:nvars-1;
    Vord=[Vord eigenVector(:,order(i,2))];    %reorder eigenvector according to ordered
end                                 %eigenvalues            
                          
%impulse vectors associated with the largest eigenvalues subject to
%zero-impact constraint
gamma = zeros(nvars,1);
gammaHat=Vord(:,1); %VHat(:,1);  % use directely Vord matrix (AK)
gamma = [0;gammaHat];   %first element of gamma = 0 so as to implement assumption that news shock has no contemp effect on TFP (ordered 1st)
A1=Gtilde*gamma;   %nvars x 1 matrix with A1 in the column 

v = gamma'*inv(Gtilde)*res';    %1xT vector of news shocks
v = v';
    
    
% Compute impulse responses
%----------------------------------------------------------------------
U1=[A1; zeros(nvars*nlags-nvars,1)];

for k=1:nhorizon;
            Zk1(k,:)=(M^(k-1)*U1)';
end
%impulse responses
impulse1(:,:)=Zk1(:,1:nvars); 
    
%in order to avert 180 degree rotations, test whether TFP impulse response is positive at some medium-run horizon (e.g. 40 quarters), otherwise revert
%ATTENTION: NEED TO CORRECTLY SPECIFY COLUMN OF RELEVANT VARIABLE; here, TFP is target for FEV maximization and TFP is 1st in VAR
if impulse1(nimp/2,1) < 0   
       A1 = - A1; 
       v = -v;
       impulse1 = -impulse1;
end


% Compute fraction of VD explained at different horizons
%-----------------------------------------------------------------------
sigmak=B(:,:,1)*vmat*B(:,:,1)';
hh1=B(:,:,1)*A1*(B(:,:,1)*A1)';
vardec1(1,:)=(diag(hh1./sigmak))';

for k=1:nhorizon-1;
        %add square of k-step ahead forecast error to build k-ahead variance-covariance 
        sigmak=sigmak+B(:,:,k+1)*vmat*B(:,:,k+1)';
        hh1=hh1+B(:,:,k+1)*A1*(B(:,:,k+1)*A1)';
        vardec1(k+1,:)=(diag(hh1./sigmak))';
end      

%simulate historical time series implied by shock
[T,nvars] = size(y);
nlags = (size(x,2)-1)/nvars;
xsim = x(1,:);        
for i=1:T
    ysim(i,:) = xsim*b + A1'*v(i);
    xsim = [ysim(i,:) xsim(1:nvars*(nlags-1)) 1];  %update vector of regressors
end
    
%build matrix of FEV and IRF outputs
data=[vardec1 impulse1];

%expand output matrix with two sheets of zeros so that we can use same plotting programs as for uhlig.m output
[temp1,temp2]=size(data);
FEVIRFs = zeros(temp1,temp2,3);  
FEVIRFs(:,:,1) = data;


