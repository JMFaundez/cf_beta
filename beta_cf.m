
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
col = linspecer(9,'sequential');

xmax = 0.3;
xmin = 0.02;
figw = 1000;
figh = 400;

figure('Position',[100 100 figw figh])
count=0;
for jj=2:10
    count = count+1;
    subplot(211)
    hold on
    yy1 = squeeze(abs(cfL_fft(jj,:)));
    yname  = "$\beta_{"+num2str(jj-1)+"}$";
    plot(squeeze(X(1,:)),yy1,'Color',col(count,:),'DisplayName',yname)
    
    subplot(212)
    hold on
    yy2 = squeeze(abs(cf_fft(jj,:)));
    yname  = "$\beta_{"+num2str(jj-1)+"}$";
    plot(squeeze(X(1,:)),yy2,'Color',col(count,:),'DisplayName',yname)
end
subplot(211)
legend("Interpreter","latex",'FontSize',14)
xlim([xmin,xmax])
xlabel('$x$','Interpreter','latex','FontSize',18)
title("Linear")

subplot(212)
xlim([xmin,xmax])
xlabel('$x$','Interpreter','latex','FontSize',18)
title("Non-Linear")


cmap = 'bluewhitered';
val = 1e-3;

figure('Position',[600 600 figw figh])
ax1 =subplot(212);
pcolor(X,Z,cf)
shading interp
clim([-val,val])
colormap(ax1,cmap)
colorbar()
axis equal
ylim([min(Z(:)),max(Z(:))])
xlim([xmin,xmax])
title('NonLinear')
ylabel('$z$','FontSize',18,'Interpreter','Latex')
xlabel('$x$','FontSize',18,'Interpreter','Latex')


ax2=subplot(211);
pcolor(X,Z,cfL)
shading interp
clim([-val,val])
colormap(ax2,cmap)
colorbar()
axis equal
ylim([min(Z(:)),max(Z(:))])
xlim([xmin,xmax])
title('Linear')
ylabel('$z$','FontSize',18,'Interpreter','Latex')





iti = 0.6/dt;
Lt = iti:2450;
Nb = 10;
Nt = length(Lt);
cf_fft = zeros(Nx,Nb,Nt);
cfL_fft = zeros(Nx,Nb,Nt);


tic
count = 0;
for it =Lt
    count = count+1;
    t = dt*it;

    [~,itNL] = min(abs(t-timeNL));
    [~,itL] = min(abs(t-timeL));

    NL = load("../NonLinear/NL_p_shear_it_"+num2str(itNL,'%5.5i'));
    L = load("../Linear/L_p_shear_it_"+num2str(itL,'%5.5i'));

    cf = (NL.cfu-BF.cfu*1)*1/Re;
    cfL = (L.cfu)*1/Re;

    cfw = (NL.cfw-BF.cfw*1)*1/Re;
    cfwL = (L.cfw)*1/Re;


    cf_ffti = fft(cf(1:end-1,:),[],1)/(Nz-1);
    cfL_ffti = fft(cfL(1:end-1,:),[],1)/(Nz-1);
    cf_fft(:,:,count) = abs(cf_ffti(2:Nb+1,:))';
    cfL_fft(:,:,count) = abs(cfL_ffti(2:Nb+1,:))';
%    for jj=1:Nb
%        cf_fft(:,jj,count) = cf_ffti(jj+1,:);
%        cfL_fft(:,jj,count) = cfL_ffti(jj+1,:);
%    end

end
toc


cf_avg = mean(cf_fft,3,'omitnan');
cfL_avg = mean(cfL_fft,3,'omitnan');

col = linspecer(Nb,'sequential');

xmax = 0.2;
xmin = 0.02;
ymin = 0;
ymax = 8e-5;
figw = 1000;
figh = 400;
lw = 1.2;

figure('Position',[100 100 figw figh])
count=0;
for jj=1:Nb
    count = count+1;
    subplot(211)
    hold on
    yy1 = squeeze(abs(cfL_avg(:,jj)));
    yname  = "$\beta_{"+num2str(jj)+"}$";
    plot(squeeze(X(1,:)),yy1,'Color',col(count,:),'LineWidth',lw,'DisplayName',yname)
    
    subplot(212)
    hold on
    yy2 = squeeze(abs(cf_avg(:,jj)));
    yname  = "$\beta_{"+num2str(jj)+"}$";
    plot(squeeze(X(1,:)),yy2,'Color',col(count,:),'LineWidth',lw,'DisplayName',yname)
end
subplot(211)
legend("Interpreter","latex",'FontSize',14)
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$c_f(\beta)$','Interpreter','latex','FontSize',18)
title("Linear")
box on

subplot(212)
xlim([xmin,xmax])
ylim([ymin,ymax])
xlabel('$x$','Interpreter','latex','FontSize',18)
ylabel('$c_f(\beta)$','Interpreter','latex','FontSize',18)
title("Non-Linear")
box on


figure('Position',[100 100 1600 800])
count=0;
for jj=1:Nb
    count = count+1;
    subplot(2,Nb/2,count)
    hold on
    yy1 = squeeze(abs(cfL_avg(:,jj)));
    yy2 = squeeze(abs(cf_avg(:,jj)));
    yname  = "$\beta_{"+num2str(jj)+"}$";
    plot(squeeze(X(1,:)),yy1,'k-','LineWidth',lw)
    plot(squeeze(X(1,:)),yy2,'r--','LineWidth',lw)
    xlim([xmin,xmax])
    ylim([ymin,ymax])
    xlabel('$x$','Interpreter','latex','FontSize',18)
    ylabel('$c_f(\beta)$','Interpreter','latex','FontSize',18)
    title(yname,'Interpreter','latex','FontSize',16)
    box on
    
end

