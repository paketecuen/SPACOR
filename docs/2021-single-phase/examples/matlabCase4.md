# Matlab simulation for case 4
## Mixed series and parallel circuit with sinusoidal source
<img src="https://electrica.ual.es/spacor/images/spacorcasefour_simulation.png" width="100%">

Matlab Code
Herein, you can find the code for simulation purposes.

```matlab                   
%**************************************************************************
%  SPACOR Space Vectors for Parameter Identification of Non-Linear Circuits
%  https://electrica.ual.es/spacor
%
%  Case 4: Mixed series and parallel circuit with sinusoidal source
%  For more information: https://electrica.ual.es/spacor/cases/index/4
%
%  Developed in MATLAB R2018a(9.4)
%
%  © SPACOR
%     @author: Francisco G. Montoya
%     @author: Francisco de León
%     @author: Francisco Arrabal-Campos 
%     @author: Alfredo Alcayde
%
%  Main paper:
%  DOI: 
%**************************************************************************

clear all;
DT = 0.0000001;
t = [0:DT:0.02-DT];
R = 0.456;
G = 1/R;
L = 1.73e-2;
Rs = 0.1;
V1 = 120;
w = 2*pi*50;
v = V1*sin(w*t);
dv = w*V1*cos(w*t);
ddv = -w^2*V1*sin(w*t);
dddv = -w^3*V1*cos(w*t);
iL = V1*w*L/(w^2*L^2+Rs^2)*exp(-Rs/L*t)-V1*w*L/(w^2*L^2+Rs^2)*cos(w*t)+V1*Rs/...
     (w^2*L^2+Rs^2)*sin(w*t);
diL = -Rs/L*V1*w*L/(w^2*L^2+Rs^2)*exp(-Rs/L*t)+V1*w^2*L/(w^2*L^2+Rs^2)*sin(w*t)...
      +V1*Rs*w/(w^2*L^2+Rs^2)*cos(w*t);
ddiL = (-Rs/L)^2*V1*w*L/(w^2*L^2+Rs^2)*exp(-Rs/L*t)+V1*w^3*L/(w^2*L^2+Rs^2)*cos(w*t)...
     -V1*Rs*w^2/(w^2*L^2+Rs^2)*sin(w*t);
dddiL = (-Rs/L)^3*V1*w*L/(w^2*L^2+Rs^2)*exp(-Rs/L*t)-V1*w^4*L/(w^2*L^2+Rs^2)*sin(w*t)...
      -V1*Rs*w^3/(w^2*L^2+Rs^2)*cos(w*t);
iR = G*v;
diR = G*dv;
ddiR = G*ddv;
dddiR = G*dddv;
ii = iR+iL;
di = diR+diL;
ddi = ddiR+ddiL;
dddi = dddiR+dddiL;
Gp = (ddi.*(v.*ddi-di.*dv)-dddi.*(v.*di-ii.*dv)+ddv.*(di.*di-ii.*ddi))...
    ./(ddi.*(v.*ddv-dv.*dv)-dddv.*(v.*di-ii.*dv)+ddv.*(dv.*di-ii.*ddv));
a  = (ddi.*(dv.*ddi-di.*ddv)-dddi.*(dv.*di-ii.*ddv)+dddv.*(di.*di-ii.*ddi))...
    ./(dddi.*(v.*ddv-dv.*dv)-dddv.*(v.*ddi-di.*dv)+ddv.*(dv.*ddi-di.*ddv));
Rser = 1./(a-Gp);
Lser = -Rser.*(ddi.*(v.*ddv-dv.*dv)-dddv.*(v.*di-ii.*dv)+ddv.*(dv.*di-ii.*ddv))...
    ./(dddi.*(v.*ddv-dv.*dv)-dddv.*(v.*ddi-di.*dv)+ddv.*(dv.*ddi-di.*ddv));

subplot(2,1,1);
hold on;
grid on;
plot(t,v,'r','LineWidth',2);
plot(t,ii,'b','LineWidth',2);
text(0.005,155,'$v(t)$','interpreter','latex','fontsize',18);
text(0.018,-220,'$i(t)$','interpreter','latex','fontsize',18);
set(gca,'FontSize',20,'TickLabelInterpreter','latex');

subplot(2,1,2);
hold on;
grid on;
plot(t,10*Rser,'r','LineWidth',2);
plot(t,500*Lser,'b','LineWidth',2);
plot(t,10./Gp,'g','LineWidth',2);
xlabel('time (s)','interpreter','latex', 'FontWeight','bold','fontsize',16);
text(0.015,0.35,'$10R_s(t)$','interpreter','latex','fontsize',18);
text(0.010,5,'$10R(t)$','interpreter','latex','fontsize',18);
text(0.002,9,'$500L(t)$','interpreter','latex','fontsize',18);
set(gca,'FontSize',20,'TickLabelInterpreter','latex');

suptitle('Case 4: Mixed series and parallel circuit with sinusoidal source');
