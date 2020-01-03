clear;close all;

load('FILD2_34559_1.900_weight_function_sm34559_1.0s.mat')
load('FILD2_34559_0-3s.mat')

pin_ngyro=FILD.pin_ngyro;
pin_npitch=FILD.pin_npitch;
scint_ngyro=FILD.scint_ngyro;
scint_npitch=FILD.scint_npitch;
pin_grid_gyro=FILD.pin_grid_gyro;
pin_grid_pitch=FILD.pin_grid_pitch;
scint_grid_gyro=FILD.scint_grid_gyro;
scint_grid_pitch=FILD.scint_grid_pitch;
W=FILD.weight_matrix_2D;

%(a,r) grid
amin = 0; %degrees
amax = 90; %degrees
rmin = 1; %cm
rmax = 6; %cm

%resolution of the tomography grid in (a,r)
adim = pin_npitch;
rdim = pin_ngyro;
da = (amax-amin)/(adim-1);
dr = (rmax-rmin)/(rdim-1);

pitch = amin:da:amax;
gyro = rmin:dr:rmax;
[R, A] = meshgrid(gyro, pitch);

%Gaussian blob in (a,r)
ac = 45;
rc = 2.5;
asigma = 1*5;
rsigma = 1;
pinhole_2D = 1 / (asigma*rsigma*2*pi) * exp(-((A-ac).^2)./(2*asigma^2) - ((R-rc).^2)./(2*rsigma^2));

pinhole_1D = reshape(pinhole_2D,pin_ngyro*pin_npitch,1);
scintillator_1D = W * pinhole_1D;

% Add normally distributed noise to scintillator_1D
N = length(scintillator_1D);
mu = mean(scintillator_1D);
sigma = std(scintillator_1D);

scintillator_1D_noise = scintillator_1D;
for i = 1:length(scintillator_1D);
    noise = normrnd(mu, sigma);
    if scintillator_1D(i) > 1e-5
        scintillator_1D_noise(i) = scintillator_1D(i) + noise;
    end
end
scintillator_2D = reshape(scintillator_1D,[scint_npitch,scint_ngyro]);
scintillator_2D_noise = reshape(scintillator_1D_noise,[scint_npitch,scint_ngyro]);

% Pinhole plot
figure(1);
imagesc(pin_grid_pitch,pin_grid_gyro,pinhole_2D');
set(gca,'ydir','normal')
colormap(hot)
set(gca, 'Fontsize', 20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
title('Pinhole distribution')
axis([30 80 1 6])

% Scintillator plot without noise
figure(2);
imagesc(scint_grid_pitch,scint_grid_gyro,scintillator_2D');
set(gca,'ydir','normal')
colormap(hot)
set(gca, 'Fontsize', 20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
title('Scintillator distribution')
axis([30 80 1 6])

% Scintillator plot with noise
figure(3);
imagesc(scint_grid_pitch,scint_grid_gyro,scintillator_2D_noise');
set(gca,'ydir','normal')
colormap(hot)
set(gca, 'Fontsize', 20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
title('Noisy scintillator distribution')
axis([30 80 1 6])

% Generate 1000 samples
samples = cell(1,1000);
targets = cell(1,1000);
noise_sigma = zeros(1,1000);

for i = 1:1000;
    i
    scintillator_1D_noise = scintillator_1D;
    sigma_val = rand();
    noise_sigma(i) = sigma_val;
    for j = 1:length(scintillator_1D);
        noise = normrnd(mu, sigma_val);
        % noise = normrnd(mu, sigma);
        if scintillator_1D(j) > 1e-5
            scintillator_1D_noise(j) = scintillator_1D(j) + noise;
        end
    end
    samples{i} = scintillator_1D_noise;
end

% Testing one of the generated samples
figure(5);
imagesc(scint_grid_pitch,scint_grid_gyro,reshape(samples{idx2},[scint_npitch,scint_ngyro])');
set(gca,'ydir','normal')
colormap(hot)
set(gca, 'Fontsize', 20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
title('Noisy scintillator distribution')
axis([30 80 1 6])


%% Tikhonov regularization and plot

alpha = linspace(0.0035, 0.35, 50);
lambda = zeros(1,1000);
for i=1:1000
    i
    [xalpha,alpha] = TikhonovNonNeg_2(W,samples{i},pin_grid_pitch,pin_grid_gyro,0);
    [d1,d2] = size(xalpha);
    normvec = zeros(1,d2);
    for j=1:d2
        normvec(j) = norm(pinhole_1D-xalpha(:,j))^2;
    end
    [val, idx] = min(normvec);
    lambda(i) = alpha(idx);
end


for i=1:d2
    normvec(i) = norm(pinhole_1D-xalpha(:,i))^2;
end
semilogx(alpha,normvec)
xlabel('\lambda')
ylabel('||F^P_{true} - F^P(\lambda)||^2')
title('Determining correct choice of \lambda')

[val_min, idx_min] = min(noise_sigma);
[val_max, idx_max] = max(noise_sigma);
% Plot for the poster - scintillator with extreme noise
figure(20);
imagesc(scint_grid_pitch,scint_grid_gyro,reshape(samples{idx_max},length(scint_grid_pitch),length(scint_grid_gyro))')
set(gca,'ydir','normal')
colormap(hot)
set(gca,'Fontsize',20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
axis([30 80 1 6])
title('Scintillator - Extreme noise level')

% Plot for the poster - scintillator with close to no noise
figure(21);
imagesc(scint_grid_pitch,scint_grid_gyro,reshape(samples{idx_min},length(scint_grid_pitch),length(scint_grid_gyro))')
set(gca,'ydir','normal')
colormap(hot)
set(gca,'Fontsize',20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
axis([30 80 1 6])
title('Scintillator - Low noise level')

% Plot for the poster - parabola shape
%plot(alpha(15:50),normvec(15:50))
semilogx(alpha(15:50),normvec(15:50))
xlabel('\lambda')
ylabel('||F^P_{true} - F^P(\lambda)||^2')
title('Determining correct choice of \lambda')

% Plot for the poster - pinhole after inversion
imagesc(pin_grid_pitch,pin_grid_gyro,reshape(xalpha(:,idx),length(pin_grid_pitch),length(pin_grid_gyro))')
set(gca,'ydir','normal')
colormap(hot)
set(gca,'Fontsize',20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
axis([30 80 1 6])

[xalpha,alpha] = TikhonovNonNeg_2(W,samples{5},pin_grid_pitch,pin_grid_gyro,0);
[d1,d2] = size(xalpha);
normvec = zeros(1,d2);
for j=1:d2
    normvec(j) = norm(pinhole_1D-xalpha(:,j))^2;
end
[val, idx] = min(normvec);
lambda(5) = alpha(idx);

% Plot for the poster - lambda as a function of noise
figure(29);clf
noise_sigma_plot = noise_sigma/sigma;
plot(noise_sigma_plot, lambda, 's')
xlabel('\sigma/\sigma_{nn} [-]')
ylabel('\lambda [-]')


%% Generate more sensible data - explore all of parameter space

% Plot - scintillator with sigma = 0.1210
figure(27);clf
imagesc(scint_grid_pitch,scint_grid_gyro,reshape(samples{915},length(scint_grid_pitch),length(scint_grid_gyro))')
set(gca,'ydir','normal')
colormap(hot)
set(gca,'Fontsize',20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
axis([30 80 1 6])
title('Scintillator \sigma = 0.1210')

% Generate 9 Gaussian blobs
ac = 35;
rc = 2.5;
asigma = 1*2;
rsigma = 1*0.2;
pinhole_2D = 1 / (asigma*rsigma*2*pi) * exp(-((A-ac).^2)./(2*asigma^2) - ((R-rc).^2)./(2*rsigma^2));

ac2 = 35;
rc2 = 4;
pinhole_2D = pinhole_2D + 1 / (asigma*rsigma*2*pi) * exp(-((A-ac2).^2)./(2*asigma^2) - ((R-rc2).^2)./(2*rsigma^2));

ac3 = 35;
rc3 = 5.5;
pinhole_2D = pinhole_2D + 1 / (asigma*rsigma*2*pi) * exp(-((A-ac3).^2)./(2*asigma^2) - ((R-rc3).^2)./(2*rsigma^2));

ac = 55;
rc = 2.5;
asigma = 1*2;
rsigma = 1*0.2;
pinhole_2D = pinhole_2D + 1 / (asigma*rsigma*2*pi) * exp(-((A-ac).^2)./(2*asigma^2) - ((R-rc).^2)./(2*rsigma^2));

ac2 = 55;
rc2 = 4;
pinhole_2D = pinhole_2D + 1 / (asigma*rsigma*2*pi) * exp(-((A-ac2).^2)./(2*asigma^2) - ((R-rc2).^2)./(2*rsigma^2));

ac3 = 55;
rc3 = 5.5;
pinhole_2D = pinhole_2D + 1 / (asigma*rsigma*2*pi) * exp(-((A-ac3).^2)./(2*asigma^2) - ((R-rc3).^2)./(2*rsigma^2));

ac = 75;
rc = 2.5;
asigma = 1*2;
rsigma = 1*0.2;
pinhole_2D = pinhole_2D + 1 / (asigma*rsigma*2*pi) * exp(-((A-ac).^2)./(2*asigma^2) - ((R-rc).^2)./(2*rsigma^2));

ac2 = 75;
rc2 = 4;
pinhole_2D = pinhole_2D + 1 / (asigma*rsigma*2*pi) * exp(-((A-ac2).^2)./(2*asigma^2) - ((R-rc2).^2)./(2*rsigma^2));

ac3 = 75;
rc3 = 5.5;
pinhole_2D = pinhole_2D + 1 / (asigma*rsigma*2*pi) * exp(-((A-ac3).^2)./(2*asigma^2) - ((R-rc3).^2)./(2*rsigma^2));

pinhole_1D = reshape(pinhole_2D,pin_ngyro*pin_npitch,1);
scintillator_1D = W * pinhole_1D;
scintillator_2D = reshape(scintillator_1D,[scint_npitch,scint_ngyro]);

% Pinhole plot
figure(31);clf
imagesc(pin_grid_pitch,pin_grid_gyro,pinhole_2D');
set(gca,'ydir','normal')
colormap(hot)
set(gca, 'Fontsize', 20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
title('Pinhole distribution')
axis([30 80 1 6])

% Scintillator plot without noise
figure(2);clf
imagesc(scint_grid_pitch,scint_grid_gyro,scintillator_2D');
set(gca,'ydir','normal')
colormap(hot)
set(gca, 'Fontsize', 20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
title('Scintillator distribution')
axis([30 80 1 6])

% % Test Tikhonov inversion
[xalpha,alpha] = TikhonovNonNeg_2(W,scintillator_1D,pin_grid_pitch,pin_grid_gyro,0);
%imagesc(pin_grid_pitch,pin_grid_gyro,reshape(xalpha(:,idx),length(pin_grid_pitch),length(pin_grid_gyro))')
imagesc(pin_grid_pitch,pin_grid_gyro,reshape(xalpha,length(pin_grid_pitch),length(pin_grid_gyro))')
set(gca,'ydir','normal')
colormap(hot)
set(gca,'Fontsize',20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
axis([30 80 1 6])

% Generate 1000 samples with noise
N_2 = length(scintillator_1D);
mu_2 = mean(scintillator_1D);
sigma_2 = std(scintillator_1D);
samples_2 = cell(1,1000);
targets_2 = cell(1,1000);
noise_sigma_2 = zeros(1,1000);

for i = 1:1000
    i
    scintillator_1D_noise = scintillator_1D;
    sigma_val = sigma_2*rand();
    noise_sigma_2(i) = sigma_val;
    for j = 1:length(scintillator_1D);
        noise = normrnd(2*mu_2, sigma_val);
        if scintillator_1D(j) > 1e-5
            scintillator_1D_noise(j) = scintillator_1D(j) + noise;
        end
    end
    samples_2{i} = scintillator_1D_noise;
end

% Testing one of the generated samples
figure(5);clf
imagesc(scint_grid_pitch,scint_grid_gyro,reshape(samples_2{983},[scint_npitch,scint_ngyro])');
set(gca,'ydir','normal')
colormap(hot)
set(gca, 'Fontsize', 20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
title('Noisy scintillator distribution')
axis([30 80 1 6])


%% Tikhonov regularization and plot

alpha = linspace(0.075, 0.13, 50);
lambda_2 = zeros(1,1000);
for i=1:1000
    i
    [xalpha,alpha] = TikhonovNonNeg_2(W,samples_2{i},pin_grid_pitch,pin_grid_gyro,0);
    [d1,d2] = size(xalpha);
    normvec = zeros(1,d2);
    for j=1:d2
        normvec(j) = norm(pinhole_1D-xalpha(:,j))^2;
    end
    [val, idx] = min(normvec);
    lambda_2(i) = alpha(idx);
end

figure(39);clf
noise_sigma_plot_2 = noise_sigma_2/sigma_2;
plot(noise_sigma_plot_2, lambda_2, 's')
xlabel('\sigma/\sigma_{nn} [-]')
ylabel('\lambda [-]')

% Exponential fit for noise_sigma_plot_2 and lambda_2
f = fit(noise_sigma_plot_2',lambda_2','poly2');
figure(40);clf;hold on
% plot(f,noise_sigma_plot_2',lambda_2')
plot(noise_sigma_plot_2',lambda_2','s')
true_lambda_2 = 0.02844*noise_sigma_plot_2.^2 + 0.002765*noise_sigma_plot_2 + 0.08266;
plot(noise_sigma_plot_2',true_lambda_2', 's')
xlabel('\sigma/\sigma_{nn} [-]')
ylabel('\lambda [-]')

% Variable with discrete lambda increments
lambda_2_disc = true_lambda_2;
idx = find(true_lambda_2 < 0.085);
lambda_2_disc(idx) = mean(true_lambda_2(idx));
idx = find(true_lambda_2 >= 0.085 & true_lambda_2 < 0.09);
lambda_2_disc(idx) = mean(true_lambda_2(idx));
idx = find(true_lambda_2 >= 0.09 & true_lambda_2 < 0.095);
lambda_2_disc(idx) = mean(true_lambda_2(idx));
idx = find(true_lambda_2 >= 0.095 & true_lambda_2 < 0.1);
lambda_2_disc(idx) = mean(true_lambda_2(idx));
idx = find(true_lambda_2 >= 0.1 & true_lambda_2 < 0.105);
lambda_2_disc(idx) = mean(true_lambda_2(idx));
idx = find(true_lambda_2 >= 0.105 & true_lambda_2 < 0.11);
lambda_2_disc(idx) = mean(true_lambda_2(idx));
idx = find(true_lambda_2 >= 0.11 & true_lambda_2 < 0.115);
lambda_2_disc(idx) = mean(true_lambda_2(idx));
figure(42);clf;hold on
plot(noise_sigma_plot_2',lambda_2_disc', 's')
xlabel('\sigma/\sigma_{nn} [-]')
ylabel('\lambda [-]')

% Test Tikhonov inversion
[xalpha,alpha] = TikhonovNonNeg_2(W,samples_2{983},pin_grid_pitch,pin_grid_gyro,0);
%imagesc(pin_grid_pitch,pin_grid_gyro,reshape(xalpha(:,idx),length(pin_grid_pitch),length(pin_grid_gyro))')
imagesc(pin_grid_pitch,pin_grid_gyro,reshape(xalpha,length(pin_grid_pitch),length(pin_grid_gyro))')
set(gca,'ydir','normal')
colormap(hot)
set(gca,'Fontsize',20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
axis([30 80 1 6])


%% Gerating new samples: 2 Gaussian blobs, 1000 measurements
% Generate 2 Gaussian blobs
ac = 35;
rc = 2.5;
asigma = 1*3;
rsigma = 1*0.4;
pinhole_2D = 1 / (asigma*rsigma*2*pi) * exp(-((A-ac).^2)./(2*asigma^2) - ((R-rc).^2)./(2*rsigma^2));

ac2 = 65;
rc2 = 4;
pinhole_2D = pinhole_2D + 1 / (asigma*rsigma*2*pi) * exp(-((A-ac2).^2)./(2*asigma^2) - ((R-rc2).^2)./(2*rsigma^2));

pinhole_1D = reshape(pinhole_2D,pin_ngyro*pin_npitch,1);
scintillator_1D = W * pinhole_1D;
scintillator_2D = reshape(scintillator_1D,[scint_npitch,scint_ngyro]);

% Pinhole plot
figure(31);clf
imagesc(pin_grid_pitch,pin_grid_gyro,pinhole_2D');
set(gca,'ydir','normal')
colormap(hot)
set(gca, 'Fontsize', 20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
title('Pinhole distribution')
axis([30 80 1 6])

% Scintillator plot without noise
figure(2);clf
imagesc(scint_grid_pitch,scint_grid_gyro,scintillator_2D');
set(gca,'ydir','normal')
colormap(hot)
set(gca, 'Fontsize', 20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
title('Scintillator distribution')
axis([30 80 1 6])

% Generate 1000 samples with noise
N_3 = length(scintillator_1D);
mu_3 = mean(scintillator_1D);
sigma_3 = std(scintillator_1D);
samples_3 = cell(1,1000);
targets_3 = cell(1,1000);
noise_sigma_3 = zeros(1,1000);

for i = 1:1000
    i
    scintillator_1D_noise_3 = scintillator_1D;
    sigma_val = 0.2*sigma_3*rand();
    noise_sigma_3(i) = sigma_val;
    for j = 1:length(scintillator_1D);
        noise = normrnd(mu_3, sigma_val);
        if scintillator_1D(j) > 0
            scintillator_1D_noise_3(j) = scintillator_1D(j) + noise;
        end
    end
    samples_3{i} = scintillator_1D_noise_3;
end

% Testing one of the generated samples
figure(5);clf
imagesc(scint_grid_pitch,scint_grid_gyro,reshape(samples_3{5},[scint_npitch,scint_ngyro])');
set(gca,'ydir','normal')
colormap(hot)
set(gca, 'Fontsize', 20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
title('Noisy scintillator distribution')
axis([30 80 1 6])


% Tikhonov regularization and plot
alpha = linspace(0.06, 0.09, 50);
lambda_3 = zeros(1,1000);
for i=1:1000
    i
    [xalpha,alpha] = TikhonovNonNeg_2(W,samples_3{i},pin_grid_pitch,pin_grid_gyro,0);
    [d1,d2] = size(xalpha);
    normvec = zeros(1,d2);
    for j=1:d2
        normvec(j) = norm(pinhole_1D-xalpha(:,j))^2;
    end
    [val, idx] = min(normvec);
    lambda_3(i) = alpha(idx);
end

figure(39);clf
noise_sigma_plot_3 = noise_sigma_3/sigma_3;
plot(noise_sigma_3, lambda_3, 's')
xlabel('\sigma/\sigma_{nn} [-]')
ylabel('\lambda [-]')

figure(40);clf
noise_sigma_plot_3 = noise_sigma_3/sigma_3;
%plot(noise_sigma_plot_3(1:10), lambda_3(1:10), 's')
plot(lambda_3, noise_sigma_3, 's')
xlabel('\sigma/\sigma_{nn} [-]')
ylabel('\lambda [-]')


% Test Tikhonov inversion
[xalpha,alpha] = TikhonovNonNeg_2(W,samples_3{1},pin_grid_pitch,pin_grid_gyro,0);
%imagesc(pin_grid_pitch,pin_grid_gyro,reshape(xalpha(:,idx),length(pin_grid_pitch),length(pin_grid_gyro))')
imagesc(pin_grid_pitch,pin_grid_gyro,reshape(xalpha,length(pin_grid_pitch),length(pin_grid_gyro))')
set(gca,'ydir','normal')
colormap(hot)
set(gca,'Fontsize',20)
xlabel('Pitch angle [^o]')
ylabel('r_L [cm]')
axis([30 80 1 6])


%% Data from ASDEX upgrade: FILD2_34559_0-3s
% This is the data that needs to be tested - no true lambdas known.
% 1. Find the 150 lambdas using the network trained
% 2. Create the movie in MATLAB using the 150 lambdas
% 3. Compare to previous movie

asdex_data = zeros(7371,150);
movieframes=1:1:150;
for itime=movieframes
    itime
    scintillator_2D=squeeze(FILD2_34559_0s_3s.data(itime,:,:))';
%    scintillator_1D=reshape(scintillator_2D',scint_ngyro*scint_npitch,1);
%    [xalpha,alpha] = TikhonovNonNeg_2(W,scintillator_1D,pin_grid_pitch,pin_grid_gyro,0);
%    [d1,d2] = size(xalpha);
%    normvec = zeros(1,d2);
%    for j=1:d2
%        normvec(j) = norm(pinhole_1D-xalpha(:,j))^2;
%    end
%    [val, idx] = min(normvec);
%    asdex_lambda(i) = alpha(idx);
    asdex_data(:,itime) = reshape(scintillator_2D',scint_ngyro*scint_npitch,1);
end

%% Creating the movie. Changed in TikhonovNonNeg_2.m function: Alpha now 
% depends on vect(itime). So ensure that vect (the y_preds from python) are
% loaded and that the running variable is called itime.

movieframes=1:1:150;
for itime=movieframes
    itime
scintillator_2D=squeeze(FILD2_34559_0s_3s.data(itime,:,:))';
scintillator_1D=reshape(scintillator_2D',scint_ngyro*scint_npitch,1);
[xalpha,alpha] = TikhonovNonNeg(W,scintillator_1D,pin_grid_pitch,pin_grid_gyro,0);
spielbergpin(:,itime)=xalpha;
spielbergscint(:,itime)=scintillator_1D;
end
save movie_asdex.mat spielbergpin spielbergscint

readNBI;
readMPl;
readMPu;

movietimes=FILD2_34559_0s_3s.time(movieframes);
movieframecount=length(movieframes);
v = VideoWriter('movie_asdex.avi');
v.FrameRate=10;
open(v);
for iitimes=11:movieframecount-25
    itime=movieframes(iitimes);
    figure(77);clf; hold on; box on;
    positionVector1 = [0.2, 0.87, 0.7, 0.1];
    subplot('Position',positionVector1)
    box on;hold on;
    %plot(MPIBu6times, MPIBu6current/11e1,'b','LineWidth',1)
    plot(NBItimes, NBIpower/1e6,'k','LineWidth',2)
    plot(movietimes(1:end), movietimes(1),'k','LineWidth',1)
    plot([movietimes(itime) movietimes(itime)], [0 20],'k','LineWidth',2)
    set(gca,'Fontsize',20)
    xlabel('time [s]')
    %ylabel('n_f [10^{18}/m^3]')
    %xlim([2.09,2.35])
    axis([movietimes(11) movietimes(end-25) 0 10])
%     set(gca,'XTick',[1.5 2.0 2.5 3.0])
%     set(gca,'XTickLabel',['1.5s'; '2.0s'; '2.5s'; '3.0s'])
    positionVector2 = [0.16, 0.15, 0.36, 0.46];
    subplot('Position',positionVector2)
    [dummy,H]=contourf(scint_grid_pitch,scint_grid_gyro,reshape(spielbergscint(:,itime),length(scint_grid_pitch),length(scint_grid_gyro))',50);
    for i = 1:length(H) set(H(i),'EdgeColor','None');end
    grid on
    grid minor
    ax=gca;
    ax.GridColor='g';
    ax.MinorGridColor='g';
    colormap(hot)
    set(gca,'Fontsize',20)
    xlabel('Pitch angle [^o]')
    ylabel('r_L [cm]')
    axis([20 80 1 6])
    caxis([0 max(max(spielbergscint))*0.15])
    title('Scintillator')
    positionVector3 = [0.59, 0.15, 0.36, 0.46];
    subplot('Position',positionVector3)
   [dummy,H]=contourf(pin_grid_pitch,pin_grid_gyro,reshape(spielbergpin(:,itime),length(pin_grid_pitch),length(pin_grid_gyro))',50);
    for i = 1:length(H) set(H(i),'EdgeColor','None');end
    grid on
    grid minor
    ax=gca;
    ax.GridColor='g';
    ax.MinorGridColor='g';
      
    colormap(hot)
    set(gca,'Fontsize',20)
    xlabel('Pitch angle [^o]')
    %ylabel('r_L [cm]')
    
    axis([20 80 1 6])
    caxis([0 max(max(spielbergscint))*0.015])
    title('Pinhole')
    writeVideo(v,getframe(gcf))
end
close(v);