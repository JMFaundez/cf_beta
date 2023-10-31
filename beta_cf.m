
clear all
addpath /scratch/josfa/Tools/matlab-tools
BF = load("../BF/BF_p_shear.mat");
mesh = load("../mesh/mesh_cf_XZ");
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

cfw = (NL.cfw-BF.cfw*1)*1/Re;
cfwL = (L.cfw)*1/Re;


cf_fft = fft(cf,[],2);


