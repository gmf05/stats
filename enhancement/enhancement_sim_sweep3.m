
% rng(0);
Ntime = 20e3;
% Nsamples = 100;
Nsamples = 100;

b0x = log(0.02);  %
b0y = log(0.01); %
bint = -1*(exp(-[1:Ntime]/50)); % intrinsic effects curve
bext = 1*(exp(-[1:Ntime]/50)); % extrinsic effects curve
xlambda = zeros(Ntime,1);
ylambda = zeros(Ntime,1);
x = zeros(Ntime, 1);
y = zeros(Ntime, 1);
xxinf = zeros(Ntime,1);
yyinf = zeros(Ntime,1);
xyinf = zeros(Ntime,1);
yxinf = zeros(Ntime,1);
xrunISI = 0;
yrunISI = 0;

% b_axis = 0.5 : 0.5 : 2;
b_axis = 0:4;
Nb = length(b_axis);
b_n = zeros(3,Nsamples);
b_all = zeros(3,Nsamples,Nb);
b_hi_n = zeros(3,Nsamples);
b_hi = zeros(3,Nsamples,Nb);
b_low_n = zeros(3,Nsamples);
b_low = zeros(3,Nsamples,Nb);
W_n = zeros(3,3,Nsamples);
W_all = zeros(3,3,Nsamples,Nb);
dev_all = zeros(4,Nsamples,Nb);
KS_all = zeros(5,Nsamples,Nb);

i1=2; b1 = b_axis(i1);
i2=3; b2 = b_axis(i2);
i3=1; b3 = b_axis(i3);

for i4 = 1:Nb

  b4 = b_axis(i4);

  [b1,b2,b3,b4]

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
    y0 = y(idx);
    Xn = X(idx,:);
    yn = y(idx);
    
    m1 = pp_model();
    m1.X = Xn(:,1); m1.y = yn; m1 = m1.fit();
    m2 = pp_model();
    m2.X = Xn(:,[1, 2]); m2.y = yn; m2 = m2.fit();
    m3 = pp_model();
    m3.X = Xn(:,[1, 3]); m3.y = yn; m3 = m3.fit();
    m4 = pp_model();
    m4.X = Xn(:,:); m4.y = yn;  m4 = m4.fit();

    b_n(:,n) = m4.b;
    conf_b = 2*sqrt(diag(m4.W));          
    b_hi_n(:,n) = m4.b + conf_b;
    b_low_n(:,n) = m4.b - conf_b;
    W_n(:,:,n) = m4.W;
    KS_all(:,n,i4) = [m1.KS(1) m2.KS(1) m3.KS(1) m4.KS(1:2)]';
    dev_all(:,n,i4) = [m1.dev m2.dev m3.dev m4.dev]';

  end

  b_all(:,:,i4) = b_n;
  b_hi(:,:,i4) = b_hi_n;
  b_low(:,:,i4) = b_low_n;
  W_all(:,:,:,i4) = W_n;
end

ratioI_all = (dev_all(3,:,:)-dev_all(4,:,:))./(dev_all(1,:,:)-dev_all(2,:,:));
ratioI_all = log(ratioI_all);
ratioX_all = (dev_all(2,:,:)-dev_all(4,:,:))./(dev_all(1,:,:)-dev_all(3,:,:));
ratioX_all = log(ratioX_all);

%save(['~/Data/enhancement_sim_sweep3.mat'], 'ratioI_all', 'ratioX_all', 'b_all',...
%  'b_hi', 'b_low', 'W_all', 'KS_all', 'dev_all', 'b_axis','Ntime','Nsamples', ...
%  'bint','bext', 'Nb');
