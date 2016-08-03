%% Seizure data 

id_type = 'LAD';
enhancement_conf = [];

for szind = 1:11
  
  %szind
  [patientName,seizureName,badChannels,szStart,szEnd,spikeParams] = getPatientInfo(szind);  
  load(['~/enhancement_permutation_sz' num2str(szind) '_LAD3.mat']);
  
  % %dev_all(:,:,bad_idx)=[]; % why doesn't this work? bc of 3d matrix?
  delta_dev = [dev_all(1,:)-dev_all(2,:);
    dev_all(1,:)-dev_all(3,:);
    dev_all(1,:)-dev_all(4,:)];
  delta_dev(isinf(delta_dev))=nan;

  enhancement_all = 1 - squeeze((delta_dev(1,:)+delta_dev(2,:))./delta_dev(3,:));
  enhancement_all(abs(enhancement_all)>=1) = nan;
  
  %
  Nchan2 = find(squeeze(dev_all(1,1,:,Nwin2))==0,1,'first');
  Nchan2 = find(squeeze(dev_all(1,1,:))==0,1,'first');
  if isempty(Nchan2), Nchan2=size(enhancement_all,2); end
  [Nchan, Nchan2]
  %[szind, Nchan2]
  %enhancement_all = enhancement_all(:,1:Nchan2,Nwin2);
  %
  
  enhancement_conf = [enhancement_conf prctile(enhancement_all, [50, 95, 5], 1)]; 
       
end

%% NJ data

load ~/Code/git/nj/NJ_info.mat;
condn = 'noise';
% condn = 'laser';

% enhancement_conf = [];
params_all = [];

% for i = 1:10
for i = 1:10
  date = info(i).date;
  
  for block = info(i).blocks
    try
    
      load(['~/Code/git/nj/enhancement/Models_' date '_B' num2str(block) '_' condn '2.mat']);
      load(['~/Data/enhancement/enhancement_permutation_' date '_B' num2str(block) '_' condn '_1000shuffles2.mat']);
      pFull = p;
      pNull = p.get_covar([1]);
      pInt = p.get_covar([1 2]);
      pExt = p.get_covar([1 3:4]);
      Nmodels = length(msNull);
      
      Ms = [msFull{:}];
      ks = [Ms.KS];
      KS = ks(1:3:end);
      bad_idx = find(KS>0.5);
            
      delta_dev = [dev_all(1,:,:)-dev_all(2,:,:);
        dev_all(1,:,:)-dev_all(3,:,:);
        dev_all(1,:,:)-dev_all(4,:,:)];
      
      enhancement_all = 1 - squeeze((delta_dev(1,:,:)+delta_dev(2,:,:))./delta_dev(3,:,:));
      %enhancement_all(:, bad_idx) = [];
      %enhancement_all(abs(enhancement_all)>=1) = nan;
      E = median(enhancement_all,1);
      
      params_all = [params_all; i*ones(Nmodels,1), ...
        block*ones(Nmodels,1), (1:Nmodels)', KS' E'];
    end
  end
end

%%
% which ones show enhancement?
enhancement_cases = find(params_all(:,5)>0);
redundant_cases = find(params_all(:,5)<-0.25);

% for i = enhancement_cases'
for i = redundant_cases(1:end)'
%   [i,find(enhancement_cases==i)]
  [i,find(redundant_cases==i)];
  
  params = params_all(i,1:3);
  date = info(params(1)).date;
  %block = info(params(1)).blocks(params(2));
  block = params(2);
  response = params(3);
  
  load(['~/Code/git/nj/enhancement/Models_' date '_B' num2str(block) '_' condn '2.mat']);
  
  pFull = p;
  pNull = p.get_covar([1]);
  pInt = p.get_covar([1 2]);
  pExt = p.get_covar([1 3:4]);
  Nmodels = length(msNull);
  m0 = msNull{response};
  m1 = msInt{response};
  m2 = msExt{response};
  m3 = msFull{response};
  
  % compare effects 
  
  figure(1);
  
  subplot(211);
  
  j=2;
  idx = pFull.covariate_ind{j}; col = 'b';
  [lags, y3] = cubic_spline(pFull.covariate_knots{j}, m3.b(idx));

  idx = pInt.covariate_ind{j}; col = 'r';
  [lags, y1] = cubic_spline(pInt.covariate_knots{j}, m1.b(idx));

  lags_ms = lags*dt*1e3;
  plot(lags_ms, exp(y3), 'b');
  hold on;
  plot(lags_ms, exp(y1), 'r');
  xlabel('Lag Time [ms]');
  ylabel('Modulation');
  title('Intrinsic')
  %xlim([0,500]);


  for j = 3
    subplot(2,1,j-1);


    idx = pFull.covariate_ind{j}; col = 'b';
    [lags, y3] = cubic_spline(pFull.covariate_knots{j}, m3.b(idx));

    idx = pExt.covariate_ind{j-1}; col = 'r';
    [lags, y2] = cubic_spline(pExt.covariate_knots{j-1}, m2.b(idx));

    lags_ms = lags*dt*1e3;
    plot(lags_ms, exp(y3), 'b');
    hold on;
    plot(lags_ms, exp(y2), 'r');
    %xlim([0,500]);
    xlabel('Lag Time [ms]');
    ylabel('Modulation');
    title('Extrinsic')

  end
  pause; clf;

  %printpdf('~/enhancement_NJ2.pdf');
  
end



%%

dt = 8/3e4;

for szind = 1:11
  
  figure(szind);
  
  [patientName,seizureName,~,szStart,szEnd] = getPatientInfo(szind);
  load(['~/enhancement_permutation_sz' num2str(szind) '_LAD4.mat']);

  m0 = msNull{1};
  m1 = msInt{1};
  m2 = msExt{1};
  m3 = msFull{1};

  subplot(211);

  j=2;
  idx = pFull.covariate_ind{j}; col = 'b';
  [lags, y3] = cubic_spline(pFull.covariate_knots{j}, m3.b(idx));

  idx = pInt.covariate_ind{j}; col = 'r';
  [lags, y1] = cubic_spline(pInt.covariate_knots{j}, m1.b(idx));

  lags_ms = lags*dt*1e3;
  plot(lags_ms, exp(y3), 'b');
  hold on;
  plot(lags_ms, exp(y1), 'r');
  xlim([0,500]);
  xlabel('Lag Time [ms]');
  ylabel('Modulation');
  title('Intrinsic')


  for j = 3
    subplot(2,1,j-1);


    idx = pFull.covariate_ind{j}; col = 'b';
    [lags, y3] = cubic_spline(pFull.covariate_knots{j}, m3.b(idx));

    idx = pExt.covariate_ind{j-1}; col = 'r';
    [lags, y2] = cubic_spline(pExt.covariate_knots{j-1}, m2.b(idx));

    lags_ms = lags*dt*1e3;
    plot(lags_ms, exp(y3), 'b');
    hold on;
    plot(lags_ms, exp(y2), 'r');
    xlim([0,50]);
    xlabel('Lag Time [ms]');
    ylabel('Modulation');
    title('Extrinsic')

  end
  pause(1);
  printpdf(['~/enhancement_SZ2_' num2str(szind) '.pdf']);
end


%% Compare cross-correlation across sim sweep
Nlags = 100;
xc = zeros(2*Nlags+1, Nsamples);
% xc = zeros(2*Nlags+1, Nsamples, Nb, Nb, Nb, Nb);

for t=1:Ntime

  xrunISI = xrunISI + 1; yrunISI = yrunISI + 1;

  yyinf(t,:) =  bint(yrunISI); xxinf(t,:) = bint(xrunISI);
  yxinf(t,:) = bext(yrunISI); xyinf(t,:) = bext(xrunISI);

  xlambda(t,:) = exp(b0x + b1*xxinf(t,:) + b2*yxinf(t,:));
  ylambda(t,:) = exp(b0y + b3*yyinf(t,:) + b4*xyinf(t,:));
  x(t,:) = poissrnd(xlambda(t,:));
  x(t, x(t,:)>=1)=1;
  xrunISI(x(t,:)==1)=0;

  y(t,:) = poissrnd(ylambda(t,:));
  y(t, y(t,:)>=1)=1;
  yrunISI(y(t,:)==1)=0;

end

for n = 1:Nsamples
  xc(:,n) = xcov(x(:,n),y(:,n),Nlags, 'coef');
end
median_xc = squeeze(prctile(xc,50,2));
lower_xc = squeeze(prctile(xc,5,2));
upper_xc = squeeze(prctile(xc,95,2));

figure
shadedErrorBar(-Nlags:Nlags, median_xc, [(upper_xc-median_xc)'; (median_xc-lower_xc)'])

%%


filepref = 'enhancement/Models';
load NJ_info.mat;
% i=1; j=1;
% date = info(i).date;
% block = info(i).blocks(j);
date = info(TASK).date;

for block = info(TASK).blocks(1)
  
  [noiseTS,laser1TS,laser2TS,laserTS] = load_trial_TS(date,block);
  trialTS = noiseTS; condn = 'noise';
  %trialTS = laserTS; condn = 'laser';
  %disp(condn);
  Ntrials = length(trialTS);
  
  d = load_pp_nj(date,block); Fs = d.Fs;
  shiftL = 0.2288; % -> window size of 30,000 @ Fs=24,414. window size 1250 @ Fs ~1kHz
  shiftR = 1;
  T = round((shiftR+shiftL)*Fs)+1; % trial length
  
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

  msNull = cell(1, d.N_channels);
  msIntrinsic = cell(1, d.N_channels);
  msExtrinsic = cell(1, d.N_channels);
  msFull = cell(1, d.N_channels);

  for response = 1

    [TASK,block,response]
    
    % set parameters
    n1=setdiff(cellsAUD,response);
    n2=setdiff(cellsPFC,response);
    p = pp_params();
    p.response = response;
    p = p.add_covar('rate',0,T_knots,T_basis);
    p = p.add_covar('intrinsic',response,Q_knots,Q_basis);
    if ~isempty(n1), p = p.add_covar('spatial-sum1',n1,R_knots,R_basis); end
    %if ~isempty(n2), p = p.add_covar('spatial-sum2',n2,R_knots,R_basis); end

    null_ind = [p.covariate_ind{1}];
    intrinsic_ind = [p.covariate_ind{[1,2]}];
    extrinsic_ind = [p.covariate_ind{[1,3:end]}];
    full_ind = [p.covariate_ind{:}];
    Ncov = p.covariate_ind{end}(end);

    if Ntrials>0

      % make design matrix
      X = zeros(T*Ntrials,Ncov);
      y = zeros(T*Ntrials,1);
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
      %X(tind(end)+1,:) = [];
      %y(tind(end)+1:end) = [];

      % fit model
      % null
      m = pp_model([filepref '-' date '-B' num2str(block) '-c' num2str(response) '-null']);
      m.fit_method = 'glmfit'; m.link='log';
      [b,dev,stats]=glmfit(X(:, null_ind),y,'poisson','constant','off');
      m.b = b; m.W = stats.covb; m.CIF=exp(X(:, null_ind)*m.b); m.y=y;
      m=m.calcGOF(); m.CIF=[]; m.y=[];
      msNull{response} = m; %m

      % intrinsic
      m = pp_model([filepref '-' date '-B' num2str(block) '-c' num2str(response) '-intrinsic']);
      m.fit_method = 'glmfit'; m.link='log';
      [b,dev,stats]=glmfit(X(:, intrinsic_ind),y,'poisson','constant','off');
      m.b = b; m.W = stats.covb; m.CIF=exp(X(:, intrinsic_ind)*m.b); m.y=y;
      m=m.calcGOF(); m.CIF=[]; m.y=[];
      msInt{response} = m; %m

      % extrinsic
      m = pp_model([filepref '-' date '-B' num2str(block) '-c' num2str(response) '-extrinsic']);
      m.fit_method = 'glmfit'; m.link='log';
      [b,dev,stats]=glmfit(X(:, extrinsic_ind),y,'poisson','constant','off');
      m.b = b; m.W = stats.covb; m.CIF=exp(X(:, extrinsic_ind)*m.b); m.y=y;
      m=m.calcGOF(); m.CIF=[]; m.y=[];
      msExt{response} = m; %m

      % full
      m = pp_model([filepref '-' date '-B' num2str(block) '-c' num2str(response) '-full']);
      m.fit_method = 'glmfit'; m.link='log';
      [b,dev,stats]=glmfit(X(:, full_ind),y,'poisson','constant','off');
      m.b = b; m.W = stats.covb; m.CIF=exp(X(:, full_ind)*m.b); m.y=y;
      m=m.calcGOF(); m.CIF=[]; m.y=[];
      msFull{response} = m; m

      d0 = msNull{response}.dev;
      d1 = msInt{response}.dev;
      d2 = msExt{response}.dev;
      d3 = msFull{response}.dev;
      enhancement = 1-(d0-d1+d0-d2)/(d0-d3)
      
      %clear X;
      %filename = [filepref '-' date '-B' num2str(block) '-' condn '-c' num2str(response) '.mat'];
    end
  end
%   filename = [filepref '_' date '_B' num2str(block) '_' condn '2.mat'];
%   save(filename,'-v7.3','msNull', 'msInt', 'msExt', 'msFull', 'p', 'd0');
%   clear msFull msNull msInt msExt
end
