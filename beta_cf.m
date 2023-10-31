
clear all
addpath /scratch/josfa/Tools/matlab-tools
BF = load("../BF/BF_p_shear.mat");
mesh = load("../mesh/mesh_cf_XZ");
X = mesh.Xc;
Z = mesh.Zc;
[Nz,Nx] = size(X);
Re = 5.33e5;
dt = 2.5e-3;
timeL = importdata("../../linear_files.dat").data(:,1);
timeNL = importdata("../../nonlinear_files.dat").data(:,1);

figure()
hold on
plot(timeL)
plot(timeNL)



it = 1200;
t = dt*it;

%itNL = find(timeNL>=t,1,'first');
%itL = find(timeL>=t,1,'first');
[~,itNL] = min(abs(t-timeNL));
[~,itL] = min(abs(t-timeL));
% figure()
% hold on
% plot(timeL,'r')
% plot(itL,timeL(itL),'ro')
% plot(timeNL,'b')
% plot(itNL,timeNL(itNL),'bo')


NL = load("../NonLinear/NL_p_shear_it_"+num2str(itNL,'%5.5i'));
L = load("../Linear/L_p_shear_it_"+num2str(itL,'%5.5i'));

cf = (NL.cfu-BF.cfu*1)*1/Re;
cfL = (L.cfu)*1/Re;

cfw = (NL.cfw-BF.cfw*1)*1/Re;
cfwL = (L.cfw)*1/Re;


cf_fft = fft(cf(1:end-1,:),[],1)/(Nz-1);
cfL_fft = fft(cfL(1:end-1,:),[],1)/(Nz-1);


beta0 = 2*pi/(Z(end,1)-Z(1,1));
betai = [0:1:15]*beta0;

figure()
for jj=2:10
    subplot(211)
    hold on
    yy1 = squeeze(abs(cfL_fft(jj,:)));
    yname  = "$\beta_{"+num2str(jj-1)+"}$";
    plot(squeeze(X(1,:)),yy1,'DisplayName',yname)
    
    subplot(212)
    hold on
    yy2 = squeeze(abs(cf_fft(jj,:)));
    yname  = "$\beta_{"+num2str(jj-1)+"}$";
    plot(squeeze(X(1,:)),yy2,'DisplayName',yname)
end
subplot(211)
legend("Interpreter","latex")
xlim([0,0.3])
xlabel('$x$','Interpreter','latex','FontSize',12)
title("Linear")

subplot(212)
xlim([0,0.3])
xlabel('$x$','Interpreter','latex','FontSize',12)
title("Non-Linear")


val = 1e-3;

figure('Position',[500 500 1000 400])
ax1 =subplot(211);
pcolor(X,Z,cf)
shading interp
clim([-val,val])
colormap(ax1,cmap)
colorbar()
axis equal
ylim([min(Z(:)),max(Z(:))])
xlim([0.02,0.5])
title('NonLinear')
ylabel('$z$','FontSize',18,'Interpreter','Latex')


ax2=subplot(212);
pcolor(X,Z,cfL)
shading interp
clim([-val,val])
colormap(ax2,cmap)
colorbar()
axis equal
ylim([min(Z(:)),max(Z(:))])
xlim([0.02,0.5])
title('Linear')
xlabel('$x$','FontSize',18,'Interpreter','Latex')
ylabel('$z$','FontSize',18,'Interpreter','Latex')
