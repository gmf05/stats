
% rng(0);
Ntime = 20e3;
% Nsamples = 1;
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

% Nb = 7;
% b_axis = linspace(0.1, 1.9, Nb);
% b_axis = linspace(0.5,5,Nb);
% b_axis = 0.5 : 0.5 : 2;
% b_axis = 0.5 : 0.5: 2;
b_axis = -5 : 5;
b_axis_pos = 1 : 5;
% b_axis = 1:4;
Nb = length(b_axis);
Nb_pos = length(b_axis_pos);
b_all = zeros(3,Nsamples,Nb_pos,Nb,Nb_pos,Nb);
b_hi = zeros(3,Nsamples,Nb_pos,Nb,Nb_pos,Nb);
b_low = zeros(3,Nsamples,Nb_pos,Nb,Nb_pos,Nb);
W_all = zeros(3,3,Nsamples,Nb_pos,Nb,Nb_pos,Nb);
dev_all = zeros(4,Nsamples,Nb_pos,Nb,Nb_pos,Nb);
KS_all = zeros(5,Nsamples,Nb_pos,Nb,Nb_pos,Nb);

i1=1;
%for i1 = 1:Nb
for i1 = 1:Nb_pos
  b1 = b_axis_pos(i1);
  i3 = i1;
  for i2 = 1:Nb
    b2 = b_axis(i2);
    %for i3 = 1:Nb_pos
      b3 = b_axis_pos(i3);
      for i4 = 1:Nb
        b4 = b_axis(i4);
        
        [i1,i2,i3,i4]
        [b1,b2,b3,b4]
        
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

          conf_b = 2*sqrt(diag(m4.W));          
          b_all(:,n,i1,i2,i3,i4) = m4.b;
          b_hi(:,n,i1,i2,i3,i4) = m4.b + conf_b;
          b_low(:,n,i1,i2,i3,i4) = m4.b - conf_b;
          W_all(:,:,n,i1,i2,i3,i4) = m4.W;
          KS_all(:,n,i1,i2,i3,i4) = [m1.KS(1) m2.KS(1) m3.KS(1) m4.KS(1:2)]';
          dev_all(:,n,i1,i2,i3,i4) = [m1.dev m2.dev m3.dev m4.dev]';
          
        end
        
      end
    %end
    
    save(['~/Data/enhancement_sim_sweep3c.mat'], 'b_all',...
      'b_hi', 'b_low', 'W_all', 'KS_all', 'dev_all', 'b_axis','Ntime','Nsamples', ...
      'bint','bext', 'Nb', 'Nb_pos', 'b_axis_pos');

  end
end

%{
%% Plot results
load(['~/Data/enhancement_sim_sweep0.mat']);

i1=8;
i2=8;
% i3=1;
% i4=1;
% ti = reshape(median(ratioI_all(:,:,4,4,:),5), [Nb,Nb]);
% ti = reshape(median(ratioI_all(4,:,:,4,:),5), [Nb,Nb]);
% ti = reshape(median(ratioI_all(:,4,4,:,:),5), [Nb,Nb]);
% ti = reshape(median(ratioI_all(4,4,:,:,:),5), [Nb,Nb]);
% ti = reshape(median(ratioI_all(:,:,3,3,:),5), [Nb,Nb]);

% ti = reshape(median(ratioI_all(:,i1,i2,:,:),1), [Nb,Nb]);
% ti = reshape(median(ratioI_all(:,i1,:,i3,:),1), [Nb,Nb]);
ti = reshape(median(ratioI_all(:,i1,i2,:,:),1), [Nb,Nb]);
% ti = reshape(median(ratioI_all(:,:,:,i3,i4),1), [Nb,Nb]);
imagesc(ti); caxis([0,2]); colorbar;

%% Find lower/upper confidence intervals

ci_all = 0*b_all;
for i1 = 1:Nb, for i2=1:Nb, for i3=1:Nb, for i4=1:Nb
  for n = 1:100,
    ci_all(:,n,i1,i2,i3,i4) = sqrt(diag(W_all(:,:,n,i1,i2,i3,i4)));
  end
end; end; end; end
b_hi = b_all + 2*ci_all;
b_low = b_all - 2*ci_all;

%% % look at confidence intervals
b_est = b_all(:,:,i1,i2,i3,i4);
b_est_hi = b_hi(:,:,i1,i2,i3,i4);
b_est_low = b_low(:,:,i1,i2,i3,i4);

k=1;
% ci = [b_est_low(k,:); b_est(k,:); b_est_hi(k,:)] % 
exp([b_est_low(k,:); b_est(k,:); b_est_hi(k,:)]) % 

%%
M=2;
% B = reshape(ratioI_all(2,2,end-M,end-M,:), 1,[]);
% B = reshape(ratioI_all(3,2,5,5,:), 1,[]);
B = reshape(ratioI_all(3,2,5,5,:), 1,[]);
B = reshape(ratioI_all(3,2,1,1,:), 1,[]);

plot(B,'bo');

[pctiles,vals] = ecdf(B);
temp = vals(getclosest(pctiles, [0.25, 0.5, 0.75]))';
hold on
plot(50*[1 1 1], temp, 'mo', 'linewidth', 10)
% ylim([0,8])
ylabel('Enhancement Ratio [Sim, Intrinsic]')
xlabel('Simulation Number')

%}


