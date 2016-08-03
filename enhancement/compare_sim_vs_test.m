
% rng(0);
Ntime = 20e3;
% % Nsamples = 1;
Nsamples = 100;

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

b_all = zeros(3,Nsamples);
b_hi = zeros(3,Nsamples);
b_low = zeros(3,Nsamples);
W_all = zeros(3,3,Nsamples);
dev_all = zeros(4,Nsamples);
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

  b = m4.b; W = m4.W;
  conf_b = 2*sqrt(diag(W));          
  b_all(:,n) = b;
  b_hi(:,n) = b + conf_b;
  b_low(:,n) = b - conf_b;
  W_all(:,:,n) = W;

  KS_all(:,n) = [m1.KS(1) m2.KS(1) m3.KS(1) m4.KS(1:2)]';
  dev_all(:,n) = [m1.dev m2.dev m3.dev m4.dev]';
end

ratioI = (dev_all(3,:)-dev_all(4,:))./(dev_all(1,:)-dev_all(2,:));
ratioX = (dev_all(2,:)-dev_all(4,:))./(dev_all(1,:)-dev_all(3,:));

distI_sim = ratioI;
distX_sim = ratioX;
dist_sim = [ratioI; ratioX];
% distI_sim = log(ratioI);
% distX_sim = log(ratioX);

% %% Reshuffling with a single trial to get enhancement ratio distribution

X1 = [ones(Ntime,1), yyinf(:,1), xyinf(:,1)];
y1 = y(:,1);

b_all = zeros(3,Nsamples);
b_hi = zeros(3,Nsamples);
b_low = zeros(3,Nsamples);
W_all = zeros(3,3,Nsamples);
dev_all = zeros(4,Nsamples);
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

  b = m4.b; W = m4.W;
  conf_b = 2*sqrt(diag(W));          
  b_all(:,n) = b;
  b_hi(:,n) = b + conf_b;
  b_low(:,n) = b - conf_b;
  W_all(:,:,n) = W;

  KS_all(:,n) = [m1.KS(1) m2.KS(1) m3.KS(1) m4.KS(1:2)]';
  dev_all(:,n) = [m1.dev m2.dev m3.dev m4.dev]';
end

ratioI = (dev_all(3,:)-dev_all(4,:))./(dev_all(1,:)-dev_all(2,:));
ratioX = (dev_all(2,:)-dev_all(4,:))./(dev_all(1,:)-dev_all(3,:));

distI_boot = ratioI;
distX_boot = ratioX;
dist_boot = [ratioI; ratioX];
% distI_boot = log(ratioI);
% distX_boot = log(ratioX);

% figure,
% subplot(211);
% boxplot([distI_sim; distI_boot]'), %ylim([0,2]);
% hold on; plot([-1 1 3], ones(1,3), 'k--', 'linewidth',3);
% subplot(212);
% boxplot([distX_sim; distX_boot]'), %ylim([0, 2]);
% hold on; plot([-1 1 3], ones(1,3), 'k--', 'linewidth',3);
% 
% %
% figure
% subplot(211), plot(1:Nsamples, [distI_sim; distX_sim]'), %ylim([0,2]);
% subplot(212), plot(1:Nsamples, [distI_boot; distX_boot]'), %ylim([0,2]);
% corrcoef(distI_sim,distX_sim)

% 2d histogram of dist_sim vs dist_boot
D1 = dist_boot; D1(D1>1e1) = nan;
D2 = dist_sim; D2(D2>1e1) = nan;

[counts1,grid] = hist3(D1');
%XY = meshgrid(grid{1},grid{2});
[counts2] = hist3(D2', grid);
figure
imagesc(grid{1},grid{2},counts1);
figure
imagesc(grid{1},grid{2},counts2);

%% Compare distribution of deviance across all resamples
% For each model







