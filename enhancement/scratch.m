rng(0);
Ntime = 1e4;
bint = -1*(exp(-[1:Ntime]/50)); % intrinsic effects curve
bext = 1*(exp(-[1:Ntime]/50)); % extrinsic effects curve
xrunISI = 0;
yrunISI = 0;
for t=1:Ntime
    xrunISI = xrunISI + 1; yrunISI = yrunISI + 1;
    yyinf(t) =  2*bint(yrunISI); xxinf(t) = 2*bint(xrunISI);
    yxinf(t) = 4*bext(yrunISI); xyinf(t) = 2*bext(xrunISI);
    xlambda(t) = exp(log(.02)+xxinf(t)+yxinf(t));
    ylambda(t) = exp(log(.01)+yyinf(t)+xyinf(t));
    x(t) = poissrnd(xlambda(t));
    if x(t)>0, x(t)=1; xrunISI = 0; end;
 
    y(t) = poissrnd(ylambda(t));
    if y(t)>0, y(t)=1; yrunISI = 0; end;
end;
 
[b0 d0 s0] = glmfit(ones(size(x')),y','poisson','constant','off');
[b1 d1 s1] = glmfit(yyinf,y','poisson');
[b2 d2 s2] = glmfit(xyinf,y','poisson');
[b3 d3 s3] = glmfit([yyinf; xyinf]',y','poisson');

ratioI = (d2-d3)/(d0-d1)
ratioX = (d1-d3)/(d0-d2)

%% FIt models, compute enhancement for 1 block

% filepref = 'enhancement/Models';
addpath('~/Code/git/nj/');
load NJ_info.mat;
i=1; j=1;
date = info(i).date;
block = info(i).blocks(j);
% date = info(TASK).date;
% block = info(TASK).blocks(1)
[noiseTS,laser1TS,laser2TS,laserTS] = load_trial_TS(date,block);
trialTS = noiseTS; condn = 'noise';
% trialTS = laserTS; condn = 'laser';
%disp(condn);
Ntrials = length(trialTS);

d = load_pp_nj(date,block); Fs = d.Fs;
shiftL = 0.2288; % -> window size of 30,000 @ Fs=24,414. window size 1250 @ Fs ~1kHz
shiftR = 1;
T = round((shiftR+shiftL)*Fs)+1; % trial length
T0 = T*Ntrials;

% select which electrodes are AUD vs PFC
% NOTE: need to do this on each block 
cellsAUD = [];
cellsPFC = [];
for n = 1:length(d.labels)
  if strfind(d.labels{n},'AUD'), cellsAUD = [cellsAUD n];
  elseif strfind(d.labels{n},'PFC'), cellsPFC = [cellsPFC n];
  else error('Electrode not in AUD or PFC'); end
end

dt_ms = round(.001 * Fs);
T_knots = [0 1]; T_basis = 'indicator'; %NJ
% T_knots = [0:0.1:1]; T_basis = 'spline'; %NJ2
% T_knots = [0 2/47 6/47 12/47 23/47 35/47 1]; T_basis = 'spline'; %NJ3
Q_knots = [0 20 50 100] * dt_ms; Q_basis = 'spline'; Q_knots(1) = 1;
R_knots = [0 10 30]  * dt_ms; R_basis = 'spline'; 
Q = length(Q_knots); R = length(R_knots);

for response = 1:d.N_channels
  response
    
  % set parameters
  n1=setdiff(cellsAUD,response);
  n2=setdiff(cellsPFC,response);
  p = pp_params();
  p.response = response;
  p = p.add_covar('rate',0,T_knots,T_basis);
  p = p.add_covar('intrinsic',response,Q_knots,Q_basis);
  if ~isempty(n1), p = p.add_covar('spatial-sum1',n1,R_knots,R_basis); end
  if ~isempty(n2), p = p.add_covar('spatial-sum2',n2,R_knots,R_basis); end

  null_ind = [p.covariate_ind{1}];
  intrinsic_ind = [p.covariate_ind{[1,2]}];
  extrinsic_ind = [p.covariate_ind{[1,3:end]}];
  full_ind = [p.covariate_ind{:}];
  Ncov = p.covariate_ind{end}(end);
  
  % make design matrix
  X = zeros(T0,Ncov);
  y = zeros(T0,1);
  tind = 0;
  for k = 1:Ntrials
    d0 = d.sub_time_fast(trialTS(k)-shiftL,trialTS(k)+shiftR).reset_time();
    d0.t = d0.t-shiftL;
    m0=pp_model();
    m0=m0.makeX(d0,p);
    tind = tind(end) + (1:d0.T);
    X(tind,:) = m0.X;
    y(tind) = d0.dn(p.response,:);
  end
  [b0,dev0,stats0]=glmfit(X(:, null_ind),y,'poisson','constant','off');
  [b1,dev1,stats1]=glmfit(X(:, intrinsic_ind),y,'poisson','constant','off');
  [b2,dev2,stats2]=glmfit(X(:, extrinsic_ind),y,'poisson','constant','off');
  [b3,dev3,stats3]=glmfit(X(:, full_ind),y,'poisson','constant','off');

  % compute deviance ratios
  ratioI_true = (dev2 - dev3) / (dev0 - dev1);
  ratioX_true = (dev1 - dev3) / (dev0 - dev2);
  [ratioI_true, ratioX_true]
end