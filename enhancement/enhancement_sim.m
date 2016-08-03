
% rng(0);
Ntime = 20e3;
dt = 1e-3;
time = (1:Ntime)*dt;

b0x = log(0.02);  %
b0y = log(0.01); %
bint = -1*(exp(-[1:Ntime]/50)); % intrinsic effects curve
bext = 1*(exp(-[1:Ntime]/50)); % extrinsic effects curve
xlambda = zeros(Ntime, 1); ylambda = zeros(Ntime, 1);
x = zeros(Ntime, 1); y = zeros(Ntime, 1);
xxinf = zeros(Ntime, 1); yyinf = zeros(Ntime, 1);
xyinf = zeros(Ntime, 1); yxinf = zeros(Ntime, 1);
xrunISI = 0; yrunISI = 0;

% b1 = 1.0;
% b2 = 1.5;
% b3 = 0.6;
% b4 = 1.5;

b1 = -1.0;
b2 = 5.0;
b3 = -1.0;
b4 = -5.0;

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


xc = xcov(x, y, 10, 'coef');

X = [ones(Ntime,1), yyinf, xyinf];

m1 = pp_model();
m1.X = X(:,1); m1.y = y; m1 = m1.fit();
m2 = pp_model();
m2.X = X(:,[1, 2]); m2.y = y; m2 = m2.fit();
m3 = pp_model();
m3.X = X(:,[1, 3]); m3.y = y; m3 = m3.fit();
m4 = pp_model();
m4.X = X(:,:); m4.y = y;  m4 = m4.fit();

b = m4.b; W = m4.W;
conf_b = 2*sqrt(diag(W));          
b_hi = b + conf_b;
b_low = b - conf_b;
% [b_low b b_hi]

KS_all = [m1.KS(1) m2.KS(1) m3.KS(1) m4.KS(1:2)]';
dev_all = [m1.dev m2.dev m3.dev m4.dev]';

delta_dev = [dev_all(1)-dev_all(2); ...
            dev_all(1)-dev_all(3); ...
            dev_all(3)-dev_all(4); ...
            dev_all(2)-dev_all(4); ...
            dev_all(1)-dev_all(4)];

enhancement = 1 - (delta_dev(1)+delta_dev(2))./delta_dev(5)

%% Compute angles like Schey 1993

% Intrinsic
% Y2 = sqrt(delta_dev(1));
% Y21 = sqrt(delta_dev(3));
% Y12 = sqrt(delta_dev(5));

% Extrinsic
Y2 = sqrt(delta_dev(2));
Y21 = sqrt(delta_dev(4));
Y12 = sqrt(delta_dev(5));

% Y21/Y12
% Y2/Y12

phi = asin(Y21/Y12);
% phi_deg = phi/pi*180;
temp = acos(Y2/Y12);
theta = temp+phi;
S = sin(phi)^2 - cos(theta-phi)^2;
[theta, phi, S]

