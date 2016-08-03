total=0;
for szind = 1:11
  [patientName,seizureName,badChannels1,szStart,szEnd,spikeParams] = getPatientInfo(szind);
  badChannels = find(sum(d.dn,2)==0);
  %badChannels = union(badChannels1, find(sum(d.dn,2)==0));
  %goodChannels = neuroport_interior_good(badChannels1);
  goodChannels = neuroport_interior_good(badChannels);
  total = total+length(goodChannels);
end

%%
FIGSIZE_CM = 4.5;
figure('units','centimeters','PaperPosition',[1 1 FIGSIZE_CM FIGSIZE_CM]);


%% Figure 2B. Plot sim results for batch

% load(['~/Data/enhancement_sim_sweep3.mat']);
% load(['~/Data/enhancement_sim_sweep3b.mat']);
load(['~/Data/enhancement_sim_sweep3c.mat']);
% load(['~/Data/enhancement_sim_sweep0.mat']);
% load(['~/Data/enhancement_sim_sweep0b.mat']);
% load(['~/Data/enhancement_sim_sweep0c.mat']);
%  load(['~/Data/enhancement_sim_sweep0d.mat']);

% for i1 = 1:3:Nb
%   for i2 = 1:3:Nb

% for i1 = 1:size(dev_all,3)
for i1 = 2
%   for i2 = 2:size(dev_all,5)
% for i1 = 4
%   for i2 = 1
    i2 = i1;
    D = squeeze(dev_all(:,:,i1,:,i2,:));
    %D = squeeze(dev_all(:,:,:,:,i1,i2));
    
    delta_dev = [D(1,:,:,:)-D(2,:,:,:);
      D(1,:,:,:)-D(3,:,:,:);
      %D(3,:,:,:)-D(4,:,:,:);
      %D(2,:,:,:)-D(4,:,:,:);
      D(1,:,:,:)-D(4,:,:,:)];
    
    enhancement = 1 - ((delta_dev(1,:,:,:)+delta_dev(2,:,:,:)) ./ delta_dev(3,:,:,:));
    median_enhancement = squeeze(median(enhancement,2));
    %lower_enhancement = squeeze(prctile(enhancement,5,2));
    %upper_enhancement = squeeze(prctile(enhancement,95,2));
    
    %median_enhancement_mask = median_enhancement;
    %median_enhancement_mask(lower_enhancement<=0)=0;
    
    figure();
    %imagesc(b_axis, b_axis, squeeze(enhancement)); caxis([-1,1]);
    imagesc(b_axis, b_axis, median_enhancement); caxis([-1,1]);
    %imagesc(median_enhancement_mask); caxis([-1,1]);
    %imagesc(b_axis_pos, b_axis, lower_enhancement); caxis([-1,1]);
    %imagesc(upper_enhancement); caxis([-1,1]);
    axis xy;
    colorbar;
    xlabel('W4');
    ylabel('W2');
    pause; %clf;
    
%   end
end

printpdf('~/enhancement_state_space2.pdf');

%% Figure 2C. Show example of bootstrapping

% load(['~/Data/enhancement_sim_sweep3.mat']);
% load(['~/Data/enhancement_sim_sweep3b.mat']);
load(['~/Data/enhancement_sim_sweep3c.mat']);
% load(['~/Data/enhancement_sim_sweep0.mat']);
% load(['~/Data/enhancement_sim_sweep0b.mat']);
% load(['~/Data/enhancement_sim_sweep0c.mat']);
% load(['~/Data/enhancement_sim_sweep0d.mat']);
% i1 = 5; i2 = 8; i3 = 8;
% i1 = 5; i2 = 1; i3 = 2;
% i1 = 2; i2 = 8; i3 = 3;
i1 = 2; i2 = 9; i3 = 2;
% i1 = 1; i2 = 9; i3 = 1;
% i1 = 2; i2 = 2; i3 = 2;
% i1 = 3; i2 = 3; i3 = 3;
% i1 = 4; i2 = 7; i3 = 5;
  
D = squeeze(dev_all(:,:,i1,i2,i3,:));
% D = squeeze(dev_all(:,:,i1,:,i3,i2));
% D = squeeze(dev_all(:,:,i1,:,i3,8));

delta_dev = [D(1,:,:,:)-D(2,:,:,:);
  D(1,:,:,:)-D(3,:,:,:);
  %D(3,:,:,:)-D(4,:,:,:);
  %D(2,:,:,:)-D(4,:,:,:);
  D(1,:,:,:)-D(4,:,:,:)];

enhancement = 1 - ((delta_dev(1,:,:)+delta_dev(2,:,:)) ./ delta_dev(3,:,:));
median_enhancement = squeeze(median(enhancement,2));
lower_enhancement = squeeze(prctile(enhancement,5,2));  
upper_enhancement = squeeze(prctile(enhancement,95,2));

figure();
%plot(b_axis, enhancement);
shadedErrorBar(b_axis, median_enhancement, [upper_enhancement'-median_enhancement'; median_enhancement'-lower_enhancement']);
hold on;
plot(b_axis([1 floor(Nb/2) Nb]), [0 0 0], 'r--', 'linewidth', 4);
ylim([-1,1]);

xlabel('W4');
ylabel('Enhancement Score');

printpdf('~/enhancement_state_space2_1d.pdf');

%% Figure 2D. Show data bootstrapping matches param bootstrapping

% Load parametric bootstrap

% load(['~/Data/enhancement_sim_sweep3.mat']);
% load(['~/Data/enhancement_sim_sweep3b.mat']);
load(['~/Data/enhancement_sim_sweep3c.mat']);
% load(['~/Data/enhancement_sim_sweep0.mat']);
% load(['~/Data/enhancement_sim_sweep0c.mat']);
% i1 = 2; i2 = 1; i3 = 2;
i1 = 1; i2 = 9; i3 = 1;
% i1 = 2; i2 = 9; i3 = 2;
% i1 = 2; i2 = 7; i3 = 3;
% i1 = 5; i2 = 8; i3 = 8;
% i1 = 2; i2 = 2; i3 = 2;
% i1 = 3; i2 = 3; i3 = 3;
% i1 = 4; i2 = 7; i3 = 5;

% Parametric bootstrapping  
D = squeeze(dev_all(:,:,i1,i2,i3,:));
%D = squeeze(dev_all(:,:,:,:,i1,i2));

delta_dev = [D(1,:,:,:)-D(2,:,:,:);
  D(1,:,:,:)-D(3,:,:,:);
  %D(3,:,:,:)-D(4,:,:,:);
  %D(2,:,:,:)-D(4,:,:,:);
  D(1,:,:,:)-D(4,:,:,:)];

enhancement_param = 1 - squeeze(((delta_dev(1,:,:)+delta_dev(2,:,:)) ./ delta_dev(3,:,:)));
enhancement_conf_param = prctile(enhancement_param, [50, 95, 5], 1);

% Do data bootstrap
Ntime = 20e3;
dt = 1e-3;
time = (1:Ntime)*dt;

b0x = log(0.02);  %
b0y = log(0.01); %
bint = -1*(exp(-[1:Ntime]/50)); % intrinsic effects curve
bext = 1*(exp(-[1:Ntime]/50)); % extrinsic e% i1 = 1; i2 = 1; i3 = 1;ffects curve
xlambda = zeros(Ntime, 1); ylambda = zeros(Ntime, 1);
x = zeros(Ntime, 1); y = zeros(Ntime, 1);
xxinf = zeros(Ntime, 1); yyinf = zeros(Ntime, 1);
xyinf = zeros(Ntime, 1); yxinf = zeros(Ntime, 1);
xrunISI = 0; yrunISI = 0;

Nsamples = 1000;
dev_all = zeros(4,Nsamples,Nb_pos,Nb,Nb_pos,Nb);

b1 = b_axis_pos(i1);
b2 = b_axis(i2);
b3 = b_axis_pos(i3);

for i4 = 1:Nb
  
  b4 = b_axis(i4);

  % Data bootstrapping
  xrunISI = 0;
  for t=1:Ntime

    xrunISI = xrunISI + 1; yrunISI = yrunISI + 1;

    yyinf(t) =  bint(yrunISI); xxinf(t) = bint(xrunISI);
    yxinf(t) = bext(yrunISI); xyinf(t) = bext(xrunISI);

    xlambda(t) = exp(b0x + b1*xxinf(t) + b2*yxinf(t));
    ylambda(t) = exp(b0y + b3*yyinf(t) + b4*xyinf(t));
    x(t) = poissrnd(xlambda(t));
    x(t, x(t)>=1)=1;
    xrunISI(x(t)==1)=0;

    y(t) = poissrnd(ylambda(t));
    y(t, y(t)>=1)=1;
    yrunISI(y(t)==1)=0;

  end
  
  X = [ones(Ntime,1), yyinf, xyinf];
  
  for n = 1:Nsamples
    %[b1,b2,b3,b4,n]
    idx = randsample(Ntime, Ntime, true);
    X0 = X(idx, :);
    y0 = y(idx);
    
    
    m = pp_model();
    m.y = y0;
    
    m1 = m; m1.X = X0(:,1); m1 = m1.fit();
    m2 = m; m2.X = X0(:,[1, 2]); m2 = m2.fit();
    m3 = m; m3.X = X0(:,[1, 3]); m3 = m3.fit();
    m4 = m; m4.X = X0(:,:); m4 = m4.fit();
    
%     conf_b = 2*sqrt(diag(m4.W));          
%     b_all(:,n,i1,i2,i3,i4) = m4.b;
%     b_hi(:,n,i1,i2,i3,i4) = m4.b + conf_b;
%     b_low(:,n,i1,i2,i3,i4) = m4.b - conf_b;
%     W_all(:,:,n,i1,i2,i3,i4) = m4.W;
%     KS_all(:,n,i1,i2,i3,i4) = [m1.KS(1) m2.KS(1) m3.KS(1) m4.KS(1:2)]';
    
    %
    dev_all(:,n,i1,i2,i3,i4) = [m1.dev m2.dev m3.dev m4.dev]';
    %
  end
end

save('param_bootstrap3c-1000.mat','dev_all','Nsamples','b_axis','b1','b2','b3','b4');

D = squeeze(dev_all(:,:,i1,i2,i3,:));
%D = squeeze(dev_all(:,:,:,:,i1,i2));

delta_dev = [D(1,:,:,:)-D(2,:,:,:);
  D(1,:,:,:)-D(3,:,:,:);
  %D(3,:,:,:)-D(4,:,:,:);
  %D(2,:,:,:)-D(4,:,:,:);
  D(1,:,:,:)-D(4,:,:,:)];

enhancement_data = 1 - squeeze(((delta_dev(1,:,:)+delta_dev(2,:,:)) ./ delta_dev(3,:,:)));
enhancement_conf_data = prctile(enhancement_data, [50, 95, 5], 1);

Ni = size(enhancement_conf_data,2);
xx = [0:10:Ni]; xx(end)=Ni;
Nx = length(xx);


figure;
shadedErrorBar(b_axis, enhancement_conf_param(1,:), [enhancement_conf_param(3,:)-enhancement_conf_param(1,:); enhancement_conf_param(1,:)-enhancement_conf_param(2,:)], 'b');
hold on;
plot(b_axis([1, floor(Ni/2), end]), [0 0 0], 'r--', 'linewidth', 4);
ylim([-1,1]);
xlabel('W4');
ylabel('Enhancement Score');
printpdf('~/enhancement_bootstrap_compare3c-1000_params.pdf');

figure;
shadedErrorBar(b_axis, enhancement_conf_data(1,:), [enhancement_conf_data(3,:)-enhancement_conf_data(1,:); enhancement_conf_data(1,:)-enhancement_conf_data(2,:)], 'r');
hold on;
plot(b_axis([1, floor(Ni/2), end]), [0 0 0], 'r--', 'linewidth', 4);
ylim([-1,1]);
xlabel('W4');
ylabel('Enhancement Score');
printpdf('~/enhancement_bootstrap_compare3c-1000_data.pdf');

% figure;
% shadedErrorBar(b_axis, enhancement_conf_data(1,:), [enhancement_conf_data(3,:)-enhancement_conf_data(1,:); enhancement_conf_data(1,:)-enhancement_conf_data(2,:)], 'k');
% hold on;
% plot(b_axis([1, floor(Ni/2), end]), [0 0 0], 'r--', 'linewidth', 4);
% ylim([-1,1]);
% xlabel('W4');
% ylabel('Enhancement Score');
% printpdf('~/enhancement_bootstrap_compare3c-100_data.pdf');

%% NJ data

load ~/Code/git/nj/NJ_info.mat;
condn = 'noise';
% condn = 'laser';

enhancement_conf = [];
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
      enhancement_all(abs(enhancement_all)>=1) = nan;
      E = median(enhancement_all,1);  
      
      enhancement_conf = [enhancement_conf prctile(enhancement_all, [50, 95, 5], 1)]; 
      params_all = [params_all; i*ones(Nmodels,1), ...
        block*ones(Nmodels,1), (1:Nmodels)', KS' E'];
    end
  end
end

[~,idx] = sort(enhancement_conf(1,:));
params_all = params_all(idx,:);
enhancement_conf = enhancement_conf(:, idx);

%%
for k = 1:20
  i=params_all(k,1);
  date = info(i).date;
  block = params_all(k,2);
  response = params_all(k,3);
  load(['~/Code/git/nj/enhancement/Models_' date '_B' num2str(block) '_' condn '2.mat']);
  load(['~/Data/enhancement/enhancement_permutation_' date '_B' num2str(block) '_' condn '_1000shuffles2.mat']);
  pFull = p;
  pNull = p.get_covar([1]);
  pInt = p.get_covar([1 2]);
  pExt = p.get_covar([1 3:4]);
  Nmodels = length(msNull);

  msFull{response}.plot(d,p);
  pause; clf;
end

%%

Ni = size(enhancement_conf,2);
xx = [0:10:Ni]; xx(end)=Ni;
Nx = length(xx);

% figure
FIGSIZE_CM = 4.5;
figure('units','centimeters','PaperPosition',[1 1 FIGSIZE_CM FIGSIZE_CM]);

shadedErrorBar(1:Ni, enhancement_conf(1,:), [enhancement_conf(3,:)-enhancement_conf(1,:); enhancement_conf(1,:)-enhancement_conf(2,:)]);
hold on; plot(xx, zeros(1,Nx), 'r--', 'linewidth', 5);
xlim([0, Ni]); ylim([-1,1]);
xlabel('Index [Sorted]');
ylabel('Enhancement Score');

% How often are cells w enhancement same cell over blocks?
ers = params_all(250:end, :);
rs = unique(params_all(250:end,[1 3]),'rows');
for r = 1:size(rs,1)
  [rs(r,:), sum(ers(:,1)==rs(r,1) & ers(:,3)==rs(r,2))]
end

printpdf('~/enhancement_NJ1.pdf');

%
% printpdf('~/enhancement_NJ1.pdf');


%% Seizure data 

id_type = 'LAD';
enhancement_conf = [];
params_all = [];

for szind = 1
  
  %szind
  [patientName,seizureName,badChannels,szStart,szEnd,spikeParams] = getPatientInfo(szind);
  filename = ['~/Data/Results/hierarchy/' patientName '_' seizureName '_pp_glm_' lower(id_type)];
  load([filename '_null.mat']);
  load([filename '_I.mat']);
  load([filename '_S.mat']);
  load([filename '_full.mat']);
  
  % % first drop bad models!
  [Nchan, Nwin] = size(msNull);
  d = load_lad2(szind);
  goodChannels = neuroport_interior_good(find(sum(d.dn')==0));
  %[szind, Nchan,length(setdiff(neuroport_interior_good(badChannels),find(sum(d.dn')==0)))]
  
  Nwin2 = floor(Nwin/2);
  bad_idx = [];
  for i = 1:Nchan
    if msFull{i,Nwin2}.KS(1)>msFull{i,Nwin2}.KS(2)
      bad_idx = [bad_idx i];
    end
  end
  %bad_idx
  
  load(['~/enhancement_permutation_sz' num2str(szind) '_LAD2.mat']);
  
  %dev_all = dev_all(:,:,:,floor(Nwin/2)); % << FOR LAD 2 only, not sz11
  delta_dev = [dev_all(1,:,:)-dev_all(2,:,:);
    dev_all(1,:,:)-dev_all(3,:,:);
    dev_all(1,:,:)-dev_all(4,:,:)];
  delta_dev(isinf(delta_dev))=nan;


  enhancement_all = 1 - squeeze((delta_dev(1,:,:)+delta_dev(2,:,:))./delta_dev(3,:,:));
  %%enhancement_all(:, bad_idx) = [];
  enhancement_all(abs(enhancement_all)>=1) = nan;
  
  enhancement_conf = [enhancement_conf prctile(enhancement_all, [50, 95, 5], 1)]; 
  params_all = [params_all; szind*ones(Nchan,1), goodChannels];
end

%%



%% For each seizure, what are top 3 electrodes for enhancement??

for szind = 1:11
  idx = params_all(:,1)==szind;
  temp = params_all(idx,2);
  [szind, temp(end-9:end)']
  
  blah = zeros(1,96);
  blah(temp(end-5:end))=1;
  plot_neuroport(blah,[]);
  title(num2str(szind));
  pause; clf;
  
end



%% Plot
% [~,idx] = sort(enhancement_conf(1,:));
% enhancement = enhancement(:, idx);
% Ni = size(enhancement,2);

[~,idx] = sort(enhancement_conf(1,:));
enhancement_conf = enhancement_conf(:, idx);
params_all = params_all(idx,:);
Ni = size(enhancement_conf,2);

% Ni = size(enhancement_conf,2);
xx = [0:10:Ni]; xx(end)=Ni;
Nx = length(xx);

FIGSIZE_CM = 4.5;
figure('units','centimeters','PaperPosition',[1 1 FIGSIZE_CM FIGSIZE_CM]);

% plot(1:Ni, enhancement(1,:));
shadedErrorBar(1:Ni, enhancement_conf(1,:), [enhancement_conf(3,:)-enhancement_conf(1,:); enhancement_conf(1,:)-enhancement_conf(2,:)]);
hold on; plot(xx, zeros(1,Nx), 'r--', 'linewidth', 5);
xlim([0, Ni]); ylim([-1,1]);
xlabel('Index [Sorted]');
ylabel('Enhancement Score');

% printpdf('~/enhancement_SZ1.pdf');

%% Seizure data model effects comparisons

id_type = 'LAD';
szind = 1;
[patientName,seizureName,badChannels,szStart,szEnd,spikeParams] = getPatientInfo(szind);
filename = ['~/Data/Results/hierarchy/' patientName '_' seizureName '_pp_glm_' lower(id_type)];
load([filename '_null.mat']);
load([filename '_I.mat']);
load([filename '_S.mat']);
load([filename '_full.mat']);
[Nchan, Nwin] = size(msNull);

% % first drop bad models!
d = load_lad2(szind);
d0 = d.downsample(8);
%[szind, Nchan,length(setdiff(neuroport_interior_good(badChannels),find(sum(d.dn')==0)))]

%% Comparing model effects

n = floor(Nwin/2);

for i = 1:Nchan

  j=2;
  
  idx = pFull.covariate_ind{j}; col = 'b';
  [lags, y] = cubic_spline(pFull.covariate_knots{j}, msFull{i,n}.b(idx));
  
  lags_ms = lags*d0.dt*1e3;
  plot(lags_ms, exp(y), col);
  hold on;
  
  idx = pInt.covariate_ind{j}; col = 'r';
  [lags, y] = cubic_spline(pInt.covariate_knots{j}, msInt{i,n}.b(idx));
  
  lags_ms = lags*d0.dt*1e3;
  plot(lags_ms, exp(y), col);
  
  pause; clf;
  
%   j=4;
%   idx = pSpace.covariate_ind{j-1}; col = 'r';
%   [lags, y] = cubic_spline(pSpace.covariate_knots{j-1}, msSpace{i,n}.b(idx));
  

end

%% Comparing model effects for simulation
% Case i.
b1 = 1; b3 = 1;
b2 = 4; b4 = -3;

% % Case ii.
% b1 = 1; b3 = 1;
% b2 = 4; b4 = 3;
% 
% % Case iii.
% b1 = 1; b3 = 1;
% b2 = -4; b4 = -3;

% rng(0);
Ntime = 20e3;
dt = 1e-3;
time = (1:Ntime)*dt;
Nsamples = 1000;

b0x = log(0.02);  %
b0y = log(0.01); %
bint = -1*(exp(-[1:Ntime]/50)); % intrinsic effects curve
bext = 1*(exp(-[1:Ntime]/50)); % extrinsic effects curve
xlambda = zeros(Ntime,Nsamples);
ylambda = zeros(Ntime,Nsamples);
x = zeros(Ntime, Nsamples);
y = zeros(Ntime, Nsamples);
xxinf = zeros(Ntime,Nsamples);
yyinf = zeros(Ntime,Nsamples);
xyinf = zeros(Ntime,Nsamples);
yxinf = zeros(Ntime,Nsamples);
xrunISI = zeros(1,Nsamples);
yrunISI = zeros(1,Nsamples);

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

b3_est = zeros(2,Nsamples);
b4_est = zeros(2,Nsamples);
b0_est = zeros(4,Nsamples);

for n = 1:Nsamples
  n
  %[b1,b2,b3,b4,n]

  Xn = [ones(Ntime,1), yyinf(:,n), xyinf(:,n)];

  m1 = pp_model();
  m1.X = Xn(:,1); m1.y = y(:,n); m1 = m1.fit();
  m2 = pp_model();
  m2.X = Xn(:,[1, 2]); m2.y = y(:,n); m2 = m2.fit();
  m3 = pp_model();
  m3.X = Xn(:,[1, 3]); m3.y = y(:,n); m3 = m3.fit();
  m4 = pp_model();
  m4.X = Xn(:,:); m4.y = y(:,n);  m4 = m4.fit();
  
%   b1_n = [m1.b - 2*sqrt(diag(m1.W)), m1.b + 2*sqrt(diag(m1.W))];
%   b2_n = [m2.b - 2*sqrt(diag(m2.W)), m2.b + 2*sqrt(diag(m2.W))];
%   b3_n = [m3.b - 2*sqrt(diag(m3.W)), m3.b + 2*sqrt(diag(m3.W))];
%   b4_n = [m4.b - 2*sqrt(diag(m4.W)), m4.b + 2*sqrt(diag(m4.W))];
  
%   b0y
%   b1_n(1,:)
%   b2_n(1,:)
%   b3_n(1,:)
%   b4_n(1,:)
  
  b0_est(:,n) = [m1.b(1); m2.b(1); m3.b(1); m4.b(1)];
  b3_est(:,n) = [m2.b(2); m4.b(2)];
  b4_est(:,n) = [m3.b(2); m4.b(3)];

end

mn3 = mean((b3_est-b3),2)
sd3 = std((b3_est-b3)')'
mn4 = mean((b4_est-b4),2)
sd4 = std((b4_est-b4)')'
save(['~/stats_' num2str(b2) '_' num2str(b4) '.mat'], 'mn3', 'sd3', 'mn4', 'sd4');


%%

load(['~/Data/enhancement_sim_sweep3c.mat']);
