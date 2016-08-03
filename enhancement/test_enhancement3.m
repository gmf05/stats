%% Permutation test for enhancement: Seizure
% szind = 1;
szind = TASK;
id_type = 'LAD';

windowSize = 30;
downsampleFactor = 8;

[patientName,seizureName,badChannels1,szStart,szEnd] = getPatientInfo(szind);
filename = [patientName '_' seizureName '_pp_glm_lad'];
d = load_lad2(szind);
d = d.downsample(downsampleFactor);

%-- declare model parameters
dt_ms = round(.001 / d.dt);
T_knots = [0 1]; T_basis = 'indicator';
Q_knots = [0 50 100 300 600 900 1200 1500:500:2500] * dt_ms; Q_basis = 'spline'; Q_knots(1) = 1; % for big spikes
R_knots = [0 1:5:51 100 200]  * dt_ms; R_basis = 'spline'; R_knots(1) = 1; % for big spikes
Q = length(Q_knots); R = length(R_knots);

pFull = pp_params();
pFull = pFull.add_covar('rate', 0, T_knots, T_basis);
pFull = pFull.add_covar('intrinsic', 0, Q_knots, Q_basis);
pFull = pFull.add_covar('spatial1', 0, R_knots, R_basis);
pFull = pFull.add_covar('spatial2', 0, R_knots, R_basis);
pFull = pFull.add_covar('spatial3', 0, R_knots, R_basis);
pFull = pFull.add_covar('spatial4', 0, R_knots, R_basis);
Ncovar = pFull.covariate_ind{end}(end);

pNull = pFull.get_covar(1);
indNull = [pFull.covariate_ind{1}];
pInt = pFull.get_covar([1, 2]);
indInt = [pFull.covariate_ind{1:2}];
pSpace = pFull.get_covar([1, 3:6]);
indSpace = [pFull.covariate_ind{[1 3:6]}];

%-- initialize arrays for storing results
% badChannels = find(sum(d.dn,2)==0);
% badChannels = union(badChannels1, find(sum(d.dn,2)==0));
goodChannels = neuroport_interior_good(badChannels);
Nchan = length(goodChannels);

winSize = 30;
startTimes = szStart : szEnd-winSize;
endTimes = startTimes + winSize;
Nwin = length(startTimes);

Nshuffles = 100;
% dev_all = zeros(4,Nshuffles,Nchan,Nwin);
dev_all = zeros(4,Nshuffles,Nchan);
msFull = cell(1,Nchan);

%%

n = floor(Nwin/2);
% for n = 1
fprintf(['\n\n Time ' num2str(startTimes(n)) ' - ' num2str(endTimes(n)) '\n\n']);
d0 = d.sub_time_fast(startTimes(n),endTimes(n));

for i = 1:Nchan
  i
  m = pp_model();
  response = goodChannels(i);
  neighbors = neuroport_neighbors(response);

  pFull.response = response;
  pFull.covariate_channels{2} = response;
  pFull.covariate_channels{3} = neighbors(1);
  pFull.covariate_channels{4} = neighbors(2);
  pFull.covariate_channels{5} = neighbors(3);
  pFull.covariate_channels{6} = neighbors(4);
  m = m.makeX(d0, pFull);
  X = m.X;
  y = d0.dn(response,:)';

  m1 = m; m1.X = X(:,indNull); m1.y = y;
  m1 = m1.fit();

  m2 = m; m2.X = X(:,indInt); m2.y = y;
  m2 = m2.fit();

  m3 = m; m3.X = X(:,indSpace); m3.y = y;
  m3 = m3.fit();

  m4 = m; m4.X = X(:,:); m4.y = y;
  m4 = m4.fit();
  m4.X = []; m4.y = []; m4.CIF = [];
  msFull{i} = m4;

  dev_true(:,i) = [m1.dev; m2.dev; m3.dev; m4.dev];
  
  if m4.KS(1)<m4.KS(2)
    for s = 1:Nshuffles
      if mod(s,10)==0, s, end
      idx = randsample(d0.T, d0.T, true);
      X0 = X(idx, :);
      y0 = y(idx);

      m1 = m; m1.X = X0(:,indNull); m1.y = y0;
      m1 = m1.fit();

      m2 = m; m2.X = X0(:,indInt); m2.y = y0;
      m2 = m2.fit();

      m3 = m; m3.X = X0(:,indSpace); m3.y = y0;
      m3 = m3.fit();

      m4 = m; m4.X = X0(:,:); m4.y = y0;
      m4 = m4.fit();

      dev_all(:,s,n) = [m1.dev; m2.dev; m3.dev; m4.dev];
    end
  end
  
  p=pFull;
  save(['~/enhancement_permutation_sz' num2str(szind) '_' id_type '4.mat'], 'dev_all', 'dev_true', 'goodChannels', 'badChannels', 'p', 'd0');
end

