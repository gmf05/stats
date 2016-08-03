#\!/bin/bash
unset DISPLAY
#Tue Jul 12 00:08:36 EDT 2016
matlab -nodisplay -singleCompThread -r "cd ..; addpath('/usr3/graduate/gfiddy/Code/matlab'); addpath('/usr3/graduate/gfiddy/Code/git/pp_tools'); run('/usr3/graduate/gfiddy/Code/git/mgh/mghtools'); pp_tools(); global TASK NTASKS; TASK = str2num(getenv('MATLAB_TASK')); NTASKS = str2num(getenv('MATLAB_NTASKS')); run('test_enhancement'); exit"

# matlab script test_enhancement.m at time of execution
: <<'COMMENT'
% --- test_enhancement.m : Tue Jul 12 00:08:36 EDT 2016 --- 

%% Permutation test for enhancement: NJ data

% filepref = 'enhancement/Models';
addpath('~/Code/git/nj/');
load NJ_info.mat;
% i=1; j=1;
% date = info(i).date;
% block = info(i).blocks(j);
date = info(TASK).date;
Nshuffles = 1e3;

for block = info(TASK).blocks
  [noiseTS,laser1TS,laser2TS,laserTS] = load_trial_TS(date,block);
  %trialTS = noiseTS; condn = 'noise';
  trialTS = laserTS; condn = 'laser';
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
  
  dev_start = zeros(4,d.N_channels);
  dev_all = zeros(4,Nshuffles,d.N_channels);
  
%   b_all = zeros(Ncov,Nsamples,Nb);
%   b_hi = zeros(Ncov,Nsamples,Nb);
%   b_low = zeros(Ncov,Nsamples,Nb);
%   W_all = zeros(Ncov,Ncov,Nsamples,Nb);
  
  for response = 1:d.N_channels

    [TASK,block,response]
    
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
    
    if Ntrials>0

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
      [b1,dev1,stats1]=glmfit(X(:, null_ind),y,'poisson','constant','off');
      [b2,dev2,stats2]=glmfit(X(:, intrinsic_ind),y,'poisson','constant','off');
      [b3,dev3,stats3]=glmfit(X(:, extrinsic_ind),y,'poisson','constant','off');
      [b4,dev4,stats4]=glmfit(X(:, full_ind),y,'poisson','constant','off');

      % compute deviance ratios
      dev_start(:,response) = [dev1; dev2; dev3; dev4];
      
      for n = 1:Nshuffles
        if mod(n,100)==0, n, end
        idx = randsample(T0, T0, true);
        X0 = X(idx, :);
        y0 = y(idx);
        
        % fit models
        [b1,dev1,stats1]=glmfit(X0(:, null_ind),y0,'poisson','constant','off');
        [b2,dev2,stats2]=glmfit(X0(:, intrinsic_ind),y0,'poisson','constant','off');
        [b3,dev3,stats3]=glmfit(X0(:, extrinsic_ind),y0,'poisson','constant','off');
        [b4,dev4,stats4]=glmfit(X0(:, full_ind),y0,'poisson','constant','off');
        
        % compute deviance ratios
%         conf_b = 2*sqrt(diag(stats.covb));
%         b_all(:,n,response) = b4;
%         b_hi(:,n,response) = b4 + conf_b;
%         b_low(:,n,response) = b4 - conf_b;
%         W_all(:,:,n,response) = stats4.covb;
        dev_all(:,n,response) = [dev1; dev2; dev3; dev4];
      end
      
    end
  end
  
  %save(['~/Data/enhancement/enhancement_permutation_' date '_B' num2str(block) '_' condn '_1000shuffles2.mat'], 'dev_start', 'dev_all', 'p', 'd0');
  save(['~/enhancement_permutation_' date '_B' num2str(block) '_' condn '_1000shuffles2.mat'], 'dev_start', 'dev_all', 'p', 'd0');
  
end

%% end of matlab script
COMMENT

