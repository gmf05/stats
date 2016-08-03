%% Compare distribution of deviance across all resamples
% For each model


%% Do it for simulated data vs. resampling
% Seizure data vs simulation etc


% rng(0);
Ntime = 20e3;
% % Nsamples = 1;
% % Nsamples = 100;
Nsamples = 5000;

b0x = log(0.01);  %
b0y = log(0.005); %
bint = -1*(exp(-[1:Ntime]/50)); % intrinsic effects curve
bext = 1*(exp(-[1:Ntime]/50)); % extrinsic effects curve

% % Multi run
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

dev_sim_all = zeros(4,Nsamples);
KS_all = zeros(5,Nsamples);

% % Enhancement (high I, moderate X)
b1 = 1.3;
b2 = 2.0;
b3 = 1.3;
b4 = 2.0;

% No enhancement
% b1 = 3.0;
% b2 = 0.0;
% b3 = 3.0;
% b4 = 0.0;
% b3 = 2.0;
% b4 = 0.3;

% % Multi run
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


% %% Gnnerating independent samples to get enhancement values

for n = 1:Nsamples
  if mod(n,10)==0, n, end
  X = [ones(Ntime,1), yyinf(:,n), xyinf(:,n)];

  m1 = pp_model();
  m1.X = X(:,1); m1.y = y(:,n); m1 = m1.fit();
  m2 = pp_model();
  m2.X = X(:,[1, 2]); m2.y = y(:,n); m2 = m2.fit();
  m3 = pp_model();
  m3.X = X(:,[1, 3]); m3.y = y(:,n); m3 = m3.fit();
  m4 = pp_model();
  m4.X = X(:,:); m4.y = y(:,n);  m4 = m4.fit();

  KS_all(:,n) = [m1.KS(1) m2.KS(1) m3.KS(1) m4.KS(1:2)]';
  dev_sim_all(:,n) = [m1.dev m2.dev m3.dev m4.dev]';
end

save(['~/Data/compare_sim_boot2.mat'], 'dev_sim_all', 'dev_boot_all', 'b1', 'b2', 'b3', 'b4');

% %% Reshuffling with a single trial to get enhancement ratio distribution

X1 = [ones(Ntime,1), yyinf(:,1), xyinf(:,1)];
y1 = y(:,1);

b_all = zeros(3,Nsamples);
b_hi = zeros(3,Nsamples);
b_low = zeros(3,Nsamples);
W_all = zeros(3,3,Nsamples);
dev_boot_all = zeros(4,Nsamples);
KS_all = zeros(5,Nsamples);

for n = 1:Nsamples
  if mod(n,10)==0, n, end
  idx = randsample(Ntime, Ntime, true);
  Xn = X1(idx,:); % reshuffle X1
  yn = y1(idx); % reshuffle y1
  
  m1 = pp_model();
  m1.X = Xn(:,1); m1.y = yn; m1 = m1.fit();
  m2 = pp_model();
  m2.X = Xn(:,[1, 2]); m2.y = yn; m2 = m2.fit();
  m3 = pp_model();
  m3.X = Xn(:,[1, 3]); m3.y = yn; m3 = m3.fit();
  m4 = pp_model();
  m4.X = Xn(:,:); m4.y = yn;  m4 = m4.fit();

  KS_all(:,n) = [m1.KS(1) m2.KS(1) m3.KS(1) m4.KS(1:2)]';
  dev_boot_all(:,n) = [m1.dev m2.dev m3.dev m4.dev]';
end

save(['~/Data/compare_sim_boot2.mat'], 'dev_sim_all', 'dev_boot_all', 'b1', 'b2', 'b3', 'b4');

% [y1,x1]=hist(dev_sim_all(1,:),6);
% [y2,x2]=hist(dev_boot_all(1,:));
% bar(x1,y1,'b'); hold on;
% bar(x2,y2,'r');

delta_dev_boot = [(dev_boot_all(1,:)-dev_boot_all(2,:)) ./ (dev_boot_all(1,:)-dev_boot_all(4,:)); ...
  (dev_boot_all(1,:)-dev_boot_all(3,:)) ./ (dev_boot_all(1,:)-dev_boot_all(4,:)); ...
  (dev_boot_all(3,:)-dev_boot_all(4,:)) ./ (dev_boot_all(1,:)-dev_boot_all(4,:)); ...
  (dev_boot_all(2,:)-dev_boot_all(4,:)) ./ (dev_boot_all(1,:)-dev_boot_all(4,:))];

delta_dev_sim = [(dev_sim_all(1,:)-dev_sim_all(2,:)) ./ (dev_sim_all(1,:)-dev_sim_all(4,:)); ...
  (dev_sim_all(1,:)-dev_sim_all(3,:)) ./ (dev_sim_all(1,:)-dev_sim_all(4,:)); ...
  (dev_sim_all(3,:)-dev_sim_all(4,:)) ./ (dev_sim_all(1,:)-dev_sim_all(4,:)); ...
  (dev_sim_all(2,:)-dev_sim_all(4,:)) ./ (dev_sim_all(1,:)-dev_sim_all(4,:))];


ratioX = (dev_boot_all(2,:)-dev_boot_all(4,:))./(dev_boot_all(1,:)-dev_boot_all(3,:));

ratioI = (dev_sim_all(3,:)-dev_sim_all(4,:))./(dev_sim_all(1,:)-dev_sim_all(2,:));
ratioX = (dev_sim_all(2,:)-dev_sim_all(4,:))./(dev_sim_all(1,:)-dev_sim_all(3,:));

ratioI_norm= (dev_sim_all(3,:)-dev_sim_all(4,:))./(dev_sim_all(1,:)-dev_sim_all(4,:));
ratioX_norm = (dev_sim_all(2,:)-dev_sim_all(4,:))./(dev_sim_all(1,:)-dev_sim_all(4,:));



%% Do bootstrapping and generating samples converge??

