%dataset_quarterly.m
%
%reads in data and transforms to desired dataset
%
%Last modified: Andre Kurmann, January 2020

function [Y, X, vars, lev]=dataset_quarterly_bigVAR(esty1,estq1,esty2,estq2,nlags)

%load data and transform
%-----------------------
dataseries = readmatrix('data_BS_VAR.xls', ...
                        'Sheet', 1, ...
                        'Range', 'B7:AM298');

varnames =   {'  Output   '      %1) real gdp per capita
               '         Consumption       '      %2) real c per capita
               '          Investment       '      %3) real i per capita
               '            Hours          '      %4) hours per capita
               '    F-R total econ hours   '      %5) hours per capita, total econ (Francis-Ramey JMCB)					                                              
               '  F-R total econ hours adj '      %6) hours per capita adj for demographic compositional changes, total econ (Francis-Ramey,JMCB)
               '         Inflation         '      %7) GDP inflation  					                                                
               '       Inflation       '      %8) CPI inflation					                                                
               '       PCE inflation       '      %9) PCE inflation					                                                
               '      PCEle infation       '      %10) PCE inflation less food & energy					                                              
               '     Federal Funds rate    '      %11) fed funds rate
               '             TFP           '      %12) Fernald's utilization-adjusted TFP series (2024 vintage)					                                                
               '    Adjusted TFP (2013)    '      %13) Fernald's utilization-adjusted TFP series (2013 vintage)
               '    Adjusted TFP (2007)    '      %14) Fernald's utilization-adjusted TFP series (2007 vintage)
               '    Utilization (2016)     '      %15) Fernald's utiliation of capital and labor (2016 vintage, using BFFK estimates)
               '    Utilization (2013)     '      %16) Fernald's utilization of capital and labor (2013 vintage, using BFK estimates)
               '    Utilization (2007)     '      %17) Fernald's utilization of capital and labor (2007 vintage, using BFK estimates)
               '   Unadjusted TFP (2016)   '      %18) Fernald's un-adjusted TFP series (2016 vintage)
               '   Unadjusted TFP (2013)   '      %19) Fernald's un-adjusted TFP series (2013 vintage)
               '   Unadjusted TFP (2007)   '      %20) Fernald's un-adjusted TFP series (2007 vintage)
               '      Aggregate Lprod      '      %21) aggregate labor productivity (real gdp (1) - hours (4) )
               '    Fernald Lprod (2016)   '      %22) Fernald's labor productivity implied by his data for k and h+LQ (2016 vintage)
               '     Capital deepening     '      %23) Fernald's capital deepending, defined as Lprod - TFP = alpha*(K-(h+LQ)) (2016 vintage)
               '   investment price index  '      %24) inverse of rel price of investment pc/pi (from Justiniano, Primiceri, Tambalotti)				                                               
               '    Consumer confidence    '      %25) Michigan cons. conf index E5Y		
               '      BOS confidence       '      %26) Bachman, Elstner, Sims BOS confidence index               
               '     Real S&P500 index     '      %27) Shiller's s&p comp index, cpi deflated, per capita					                                                
               '  GRR return on capital    '      %28) Gomme, Ravikumar and Rupert's return on capital
               '     R&D expenditures      '      %29) Real expenditures in Research & Development
               '     5yr - ffr spread      '      %30) 5-year - ffr bond yield spread
               ' Fama-French Ex Mkt return '      %31) Fama-French Excess Market return
               'Baron-Schmidt SSO ICT count'      %32) SSO ICT baseline count (in logs)
               'Baron-Schmidt SSO ICT+ELEC '      %33) SSO ICT+ELEC count (in logs)
               ' JLN macro uncertainty u1  '      %34) Jurado, Ludvigson, Ng macro uncertainty 1-month out
               ' JLN macro uncertainty u3  '      %35) ...3 months out
               ' JLN macro uncertainty u12 '      %36) ...12 months out
               '       Vacancies           '      %37) Vacancies
               '       Unemployment        '      %38) Unemployment Level (not rate)
               
               
               
               };   

               
%rescaling variables that are in log levels instead of percent       
   dataseries=100*dataseries;
   dataseries(:,7:11)=dataseries(:,7:11)/100;
   dataseries(:,28)=dataseries(:,28)/100; 
   dataseries(:,30:31)=dataseries(:,30:31)/100;
   dataseries(:,11) = dataseries(:,11)/4; 
        
%Step 3: setup VAR and construct lags
%------------------------------------
    
var_select=[12 2 1 37 38 8  11];  %small VAR with 2016 TFPu vintage (for Figure 2 and 4)
%var_select=[14 2 4 7];  %small VAR with 2007 TFPu vintage (for Figure 2 and 4)
lev = [1 1 1 1 1 0 0];        %set Minnesota priors for first lag of each variable (1 if random walk; 0 if white noise)

% var_select=[21 2 23 7];    %small VAR with alternative productivity measures (for Figures in appendix)
% lev = [1 1 1 0];          %set Minnesota priors for first lag of each variable (1 if random walk; 0 if white noise)

 %var_select = [12 27 2 1 3 4 7 11]; %large VAR with 2016 TFPu vintage (for Figure 5 and 6)
% var_select = [14 27 2 1 3 4 7 11]; %large VAR with 2007 TFPu vintage (for Figure 5 and 6)
 %lev = [1 1 1 1 1 0 0 0];            %set Minnesota priors for first lag of each variable (1 if random walk; 0 if white noise)

% for IRFs in Figure 6, sequentially replace last variable (11) with 33, 30, 25, 31, 26, 27
      

% Set start of sample date
daty1=1947;       % First Year of Data Set
datq1=1;          % First Quarter of Data Set

%select sample period
n1 = (esty1-daty1)*4 + ( estq1 - datq1 + 1);      
n2 = (esty2-daty1)*4 + ( estq2 - datq1 + 1);

data = dataseries(n1:n2,var_select);
vars = varnames(var_select,:);

%checking for out-of-sample values
if sum(sum(isnan(data)))>0 | sum(sum(data<-9999))>0;
    data(1:5,:)
    disp('dataseries out of range')
    disp('hit enter to continue or ctrl-c to abort')
    return
end

[T,nvars]=size(data);
%compute nlags lags
for p=1:nlags;
    X(:,1+(p-1)*nvars:p*nvars)=data((nlags+1-p):(T-p),:);
    %gives (T-nlags) x nvars*nlags matrix of lags for all variables
    %first lag of all variables first, then second lag of all variables
    %and so on...                                               
end;

%resizing sample since we loose nlags observations through the lagging of
%variables for VAR part
Y=data((nlags+1):T,:);

