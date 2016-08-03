FIGSIZE_CM = 4.5;

%% Seizure data

id_type = 'SAD';
ratioI_all = [];
ratioX_all = [];
for szind = 1:11
  szind
  [patientName,seizureName,badChannels,szStart,szEnd,spikeParams] = getPatientInfo(szind);
  filename = ['~/Data/Results/hierarchy/' patientName '_' seizureName '_pp_glm_' lower(id_type)];
  load([filename '_null.mat']);
  load([filename '_I.mat']);
  load([filename '_S.mat']);
  load([filename '_full.mat']);

  [Nchan, Nwin] = size(msNull);
  dev_01 = zeros(Nchan,Nwin);
  dev_02 = zeros(Nchan,Nwin);
  dev_12 = zeros(Nchan,Nwin);
  dev_21 = zeros(Nchan,Nwin);
  for n = 1:Nchan
    for t = 1:Nwin
      try
        dev_01(n,t) = msNull{n,t}.dev - msInt{n,t}.dev;
        dev_02(n,t) = msNull{n,t}.dev - msSpace{n,t}.dev;
        dev_12(n,t) = msInt{n,t}.dev - msFull{n,t}.dev;
        dev_21(n,t) = msSpace{n,t}.dev - msFull{n,t}.dev;
      catch
        dev_01(n,t) = nan;
        dev_02(n,t) = nan;
        dev_12(n,t) = nan;
        dev_21(n,t) = nan;
      end
    end
  end
  ratioI_all = [ratioI_all; dev_21(:)./dev_01(:)];
  ratioX_all = [ratioX_all; dev_12(:)./dev_02(:)];
end

% %% Plot

figure('units','centimeters','PaperPosition',[1 1 FIGSIZE_CM FIGSIZE_CM]);
% P=ratioI_all; my_effects = 'Intrinsic'; 
P=ratioX_all; my_effects = 'Extrinsic';
P(isnan(P))=[];
P(isinf(P))=[];
P(P==0)=[];
Pmax = 2.3;
hist(P,-0.5:0.05:Pmax);
xlim([-0.1, Pmax+0.2]);
xlabel(['Enhancement Ratio [' my_effects ' Effects]'])
ylabel('Count')
print(['~/enhancement_sz' id_type '_' my_effects '.pdf'], '-dpdf');

figure('units','centimeters','PaperPosition',[1 1 FIGSIZE_CM FIGSIZE_CM]);
P=ratioI_all; my_effects = 'Intrinsic'; 
% P=ratioX_all; my_effects = 'Extrinsic';
P(isnan(P))=[];
P(isinf(P))=[];
P(P==0)=[];
Pmax = 7.5;
hist(P,-0.5:0.1:Pmax);
xlim([-0.1, Pmax+0.2]);
xlabel(['Enhancement Ratio [' my_effects ' Effects]'])
ylabel('Count')
print(['~/enhancement_sz' id_type '_' my_effects '.pdf'], '-dpdf');

%% Seizure data 2


id_type = 'LAD';
ratioI_all = [];
ratioX_all = [];
for szind = 1
  szind
  [patientName,seizureName,badChannels,szStart,szEnd,spikeParams] = getPatientInfo(szind);
  filename = ['~/Data/Results/hierarchy/' patientName '_' seizureName '_pp_glm_' lower(id_type)];
  load([filename '_null.mat']);
  load([filename '_I.mat']);
  load([filename '_S.mat']);
  load([filename '_full.mat']);

  [Nchan, Nwin] = size(msNull);
  dev_01 = zeros(Nchan,Nwin);
  dev_02 = zeros(Nchan,Nwin);
  dev_12 = zeros(Nchan,Nwin);
  dev_21 = zeros(Nchan,Nwin);
  for n = 1:Nchan
    for t = 1:Nwin
      try
        dev_01(n,t) = msNull{n,t}.dev - msInt{n,t}.dev;
        dev_02(n,t) = msNull{n,t}.dev - msSpace{n,t}.dev;
        dev_13(n,t) = msInt{n,t}.dev - msFull{n,t}.dev;
        dev_23(n,t) = msSpace{n,t}.dev - msFull{n,t}.dev;
        
        
        
      catch
        dev_01(n,t) = nan;
        dev_02(n,t) = nan;
        dev_13(n,t) = nan;
        dev_23(n,t) = nan;
      end
    end
  end
  ratioI_all = [ratioI_all; dev_23(:)./dev_01(:)];
  ratioX_all = [ratioX_all; dev_13(:)./dev_02(:)];
end

% %% Plot

[~,idx_I] = sort(ratioI_all);
all_I = ratioI_all(idx_I);
[~,idx_X] = sort(ratioX_all);
all_X = ratioX_all(idx_X);

% Drop a couple models where bounds -> Inf
cols = find(all_I==0 | all_X==0 | isinf(all_I) | isinf(all_X) | isnan(all_I) | isnan(all_X) );
all_I(cols)=[];
all_X(cols)=[];

Ni = size(all_I,1);
xx = [0:10:Ni]; xx(end)=Ni;
Nx = length(xx);

figure
plot(1:Ni, all_I);
% shadedErrorBar(1:size(all_I,2), all_I(1,:), [all_I(2,:)-all_I(1,:); all_I(1,:)-all_I(3,:)]);
hold on; plot(xx, 1*ones(1,Nx), 'r--', 'linewidth', 5);
xlim([0, Ni]); ylim([0,8]);
% xlim([-100,0]+Ni); ylim([0.7, 2]);
xlabel('Index [Sorted]');
ylabel('Enhancement Ratio [Intrinsic Effects]');

figure
plot(1:Ni, all_X);
% shadedErrorBar(1:size(all_X,2), all_X(1,:), [all_X(2,:)-all_X(1,:); all_X(1,:)-all_X(3,:)]);
hold on; plot(xx, 1*ones(1,Nx), 'r--', 'linewidth', 5);
xlim([0, Ni]); ylim([0,2]);
% xlim([-100,0]+Ni); ylim([0.9, 1.3]);
xlabel('Index [Sorted]');
ylabel('Enhancement Ratio [Extrinsic Effects]');

%% NJ data

load ~/Code/git/nj/NJ_info.mat;
condn = 'noise';
% condn = 'laser';


all_I = [];
all_X = [];
params_I = [];
params_X = [];

for i = [1:3 5:10]
  date = info(i).date;
  
  for block = info(i).blocks
    try
      load(['~/Code/git/nj/enhancement/Models_' date '_B' num2str(block) '_' condn '.mat'], 'dev_ratio_I', 'dev_ratio_X');
      load(['~/enhancement_permutation_' date '_B' num2str(block) '_' condn '_1000shuffles.mat']);
      
      [sort_I,idx_I] = sort(dev_ratio_I);
      mean_I = 0*sort_I;
      upper_I = 0*sort_I;
      lower_I = 0*sort_I;
      [sort_X,idx_X] = sort(dev_ratio_X);
      mean_X = 0*sort_X;
      upper_X = 0*sort_X;
      lower_X = 0*sort_X;
      
      for n = 1:size(ratioI,1)
        
        
        cols = find(ratioI(idx_I(n),:)>0 & ratioX(idx_X(n),:)>0);        
        %cols = find(~isinf(ratioI(idx_I(n),:)) & ~isinf(ratioX(idx_X(n),:)));
        %cols = find(ratioI(idx_I(n),:)>0 & ratioX(idx_X(n),:)>0 & ~isinf(ratioI(idx_I(n),:)) & ~isinf(ratioX(idx_X(n),:)));
        %M=size(ratioI,2)-length(cols);
        %if M>0, M, end
        
        [pctiles,vals] = ecdf(ratioI(idx_I(n),cols));
        %temp = vals(getclosest(pctiles, [0.05, 0.5, 0.95]))';
        temp = vals(getclosest(pctiles, [0.25, 0.5, 0.75]))';
        lower_I(n) = temp(1);
        mean_I(n) = temp(2);
        upper_I(n) = temp(3);
        
        [pctiles,vals] = ecdf(ratioX(idx_X(n),cols));
        %temp = vals(getclosest(pctiles, [0.05, 0.5, 0.95]))';
        temp = vals(getclosest(pctiles, [0.25, 0.5, 0.75]))';
        lower_X(n) = temp(1);
        mean_X(n) = temp(2);
        upper_X(n) = temp(3);
        
      end
      
      
      
      Nc = size(ratioI,1);
      params_I = [params_I [i*ones(1,Nc); block*ones(1,Nc); 1:Nc; mean_I; upper_I; lower_I]];
      params_X = [params_X [i*ones(1,Nc); block*ones(1,Nc); 1:Nc; mean_X; upper_X; lower_X]];
      
      all_I = [all_I [mean_I; upper_I; lower_I]];
      all_X = [all_X [mean_X; upper_X; lower_X]];
      
    end
  end
end

[~,idx_I] = sort(all_I(1,:));
all_I = all_I(:, idx_I);
[~,idx_X] = sort(all_X(1,:));
all_X = all_X(:, idx_X);

% Drop a couple models where bounds -> Inf
cols = find(all_X(2,:)>1e10);
all_I(:, cols)=[];
all_X(:, cols)=[];

Ni = size(all_I,2);
xx = [0:10:Ni]; xx(end)=Ni;
Nx = length(xx);

figure
shadedErrorBar(1:size(all_I,2), all_I(1,:), [all_I(2,:)-all_I(1,:); all_I(1,:)-all_I(3,:)]);
hold on; plot(xx, 1*ones(1,Nx), 'r--', 'linewidth', 5);
xlim([0, Ni]); ylim([0,2.5]);
% xlim([-100,0]+Ni); ylim([0.7, 2]);
xlabel('Index [Sorted]');
ylabel('Enhancement Ratio [Intrinsic Effects]');

figure
shadedErrorBar(1:size(all_X,2), all_X(1,:), [all_X(2,:)-all_X(1,:); all_X(1,:)-all_X(3,:)]);
hold on; plot(xx, 1*ones(1,Nx), 'r--', 'linewidth', 5);
xlim([0, Ni]); ylim([0,1.5]);
% xlim([-100,0]+Ni); ylim([0.9, 1.3]);
xlabel('Index [Sorted]');
ylabel('Enhancement Ratio [Extrinsic Effects]');

