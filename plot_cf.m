clear all
addpath /scratch/josfa/Tools/matlab-tools
BF = load("BF/BF_p_shear.mat");
mesh = load("mesh/mesh_cf_XZ");
X = mesh.Xc;
Z = mesh.Zc;
Re = 5.33e5;
dt = 2.5e-3;
timeL = importdata("../linear_files.dat").data(:,1);
timeNL = importdata("../nonlinear_files.dat").data(:,1);

figure()
hold on
plot(timeL)
plot(timeNL)

%% Data Diego
tic
LD = load("Linear_Diego/wallStresses.mat",'Q','t');
toc


tic
NLD = load("NonLinear_Diego/wallStresses.mat",'Q','t');
toc

%%
gridDiego = load("Linear_Diego/grid");
time_diego = LD.t;
time_diegoNL = LD.t;

[Zd,Xd] = meshgrid(gridDiego.Z,gridDiego.X(:,1)); 
[Np,Nt] = size(LD.Q);
[Ns,Nz] = size(Xd);
[~,NtNL] = size(LD.Q);

%%

it = 1200;
t = dt*it;

%itNL = find(timeNL>=t,1,'first');
%itL = find(timeL>=t,1,'first');
[~,itNL] = min(abs(t-timeNL));
[~,itL] = min(abs(t-timeL));
[~,itDL] = min(abs(t-time_diego));
[~,itDNL] = min(abs(t-time_diegoNL));
% figure()
% hold on
% plot(timeL,'r')
% plot(itL,timeL(itL),'ro')
% plot(timeNL,'b')
% plot(itNL,timeNL(itNL),'bo')


NL = load("NonLinear/NL_p_shear_it_"+num2str(itNL,'%5.5i'));
L = load("Linear/L_p_shear_it_"+num2str(itL,'%5.5i'));

cf = (NL.cfu-BF.cfu*1)*1/Re;
cfL = (L.cfu)*1/Re;

cfLD = (LD.Q(Np/3+1:2*Np/3,itDL));
cfLD = reshape(cfLD,Ns,Nz);

cfNLD = (NLD.Q(Np/3+1:2*Np/3,itDNL));
cfNLD = reshape(cfNLD,Ns,Nz);


cfw = (NL.cfw-BF.cfw*1)*1/Re;
cfwL = (L.cfw)*1/Re;

cmap = 'bluewhitered';
val = 1e-3;
valw = 1e-3;


figure('Position',[500 500 1000 400])
ax1 =subplot(211);
pcolor(Xd,Zd,cfLD)
shading interp
clim([-val,val])
colormap(ax1,cmap)
colorbar()
axis equal
ylim([min(Z(:)),max(Z(:))])
xlim([0.02,0.5])
title('Diego Linear')
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
title('Jose Linear')
xlabel('$x$','FontSize',18,'Interpreter','Latex')
ylabel('$z$','FontSize',18,'Interpreter','Latex')




figure('Position',[500 500 1000 400])
ax1 =subplot(211);
pcolor(Xd,Zd,cfNLD)
shading interp
clim([-val,val])
colormap(ax1,cmap)
colorbar()
axis equal
ylim([min(Z(:)),max(Z(:))])
xlim([0.02,0.5])
title('Diego NonLinear')
ylabel('$z$','FontSize',18,'Interpreter','Latex')


ax2=subplot(212);
pcolor(X,Z,cf)
shading interp
clim([-val,val])
colormap(ax2,cmap)
colorbar()
axis equal
ylim([min(Z(:)),max(Z(:))])
xlim([0.02,0.5])
title('Jose NonLinear')
xlabel('$x$','FontSize',18,'Interpreter','Latex')
ylabel('$z$','FontSize',18,'Interpreter','Latex')


%%

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




figure('Position',[500 500 1000 400])
ax1 =subplot(211);
pcolor(X,Z,cfw)
shading interp
clim([-valw,valw])
colormap(ax1,cmap)
colorbar()
axis equal
ylim([min(Z(:)),max(Z(:))])
xlim([0.02,0.5])
title('NonLinear')
ylabel('$z$','FontSize',18,'Interpreter','Latex')


ax2=subplot(212);
pcolor(X,Z,cfwL)
shading interp
clim([-valw,valw])
colormap(ax2,cmap)
colorbar()
axis equal
ylim([min(Z(:)),max(Z(:))])
xlim([0.02,0.5])
title('Linear')
xlabel('$x$','FontSize',18,'Interpreter','Latex')
ylabel('$z$','FontSize',18,'Interpreter','Latex')



