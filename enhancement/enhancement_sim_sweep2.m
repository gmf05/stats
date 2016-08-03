
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

b_axis = -5 : 5;
Nb = length(b_axis);
b_all = zeros(3,Nsamples,Nb,Nb);
b_hi = zeros(3,Nsamples,Nb,Nb);
b_low = zeros(3,Nsamples,Nb,Nb);
W_all = zeros(3,3,Nsamples,Nb,Nb);
dev_all = zeros(4,Nsamples,Nb,Nb);
KS_all = zeros(5,Nsamples,Nb,Nb);

b1 = -1; b3 = -1;
%for i1 = 1:Nb
for i2 = 1:Nb
  b2 = b_axis(i2);
    for i4 = 1:Nb
      b4 = b_axis(i4);
      
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
        b_all(:,n,i2,i4) = m4.b;
        b_hi(:,n,i2,i4) = m4.b + conf_b;
        b_low(:,n,i2,i4) = m4.b - conf_b;
        W_all(:,:,n,i2,i4) = m4.W;
        KS_all(:,n,i2,i4) = [m1.KS(1) m2.KS(1) m3.KS(1) m4.KS(1:2)]';
        dev_all(:,n,i2,i4) = [m1.dev m2.dev m3.dev m4.dev]';

      end

    end
    
    save(['~/Data/enhancement_sim_sweep2.mat'], 'b_all',...
      'b_hi', 'b_low', 'W_all', 'KS_all', 'dev_all', 'b_axis','Ntime','Nsamples', ...
      'bint','bext', 'Nb', 'Nb_pos', 'b_axis_pos');

end



