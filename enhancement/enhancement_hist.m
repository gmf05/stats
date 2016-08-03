FIGSIZE_CM = 4.5;

%% Seizure data

id_type = 'LAD';
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

%% NJ data

load ~/Code/git/nj/NJ_info.mat;
condn = 'noise';
% condn = 'laser';
ratioI_all = [];
ratioX_all = [];
count=0;
params_all = [];
for i = 1:10
  date = info(i).date;
  for block = info(i).blocks
    try
      load(['~/Code/git/nj/enhancement/Models_' date '_B' num2str(block) '_' condn '.mat']);
      N_channels = length(dev_ratio_I);
      %d = load_pp_nj(date,block); Fs = d.Fs;
      %ratioI_all = [ratioI_all dev_ratio_I];
      %ratioX_all = [ratioX_all dev_ratio_X];
      params_all = [params_all; i*ones(N_channels,1), block*ones(N_channels,1), (1:N_channels)' dev_ratio_I' dev_ratio_X'];
      count = count+N_channels;
    catch
      [];
    end
  end
end

%% Plot

figure('units','centimeters','PaperPosition',[1 1 FIGSIZE_CM FIGSIZE_CM]);
% P=params_all(:,4); my_effects = 'Intrinsic'; 
P=params_all(:,5); my_effects = 'Extrinsic';
P(isnan(P))=[];
P(isinf(P))=[];
P(P==0)=[];
hist(P,-0.5:0.05:2.5);
xlim([-0.1,2.2]);
xlabel(['Enhancement Ratio [' my_effects ' Effects]'])
ylabel('Count')
print(['~/enhancement_nj_' my_effects '.pdf'], '-dpdf');

figure('units','centimeters','PaperPosition',[1 1 FIGSIZE_CM FIGSIZE_CM]);
P=params_all(:,4); my_effects = 'Intrinsic'; 
% P=params_all(:,5); my_effects = 'Extrinsic';
P(isnan(P))=[];
P(isinf(P))=[];
P(P==0)=[];
hist(P,-0.5:0.05:2.5);
xlim([-0.1,2.2]);
xlabel(['Enhancement Ratio [' my_effects ' Effects]'])
ylabel('Count')
print(['~/enhancement_nj_' my_effects '.pdf'], '-dpdf');

%%


