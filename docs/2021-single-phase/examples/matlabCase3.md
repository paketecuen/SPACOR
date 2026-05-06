# Matlab simulation for case 3
## Serie RL with Sinusoidal source
<img src="https://electrica.ual.es/spacor/images/spacorcasethree_simulation.png" width="100%">

Matlab Code
Herein, you can find the code for simulation purposes.

```matlab                   
%**************************************************************************
%  SPACOR Space Vectors for Parameter Identification of Non-Linear Circuits
%  https://electrica.ual.es/spacor
%
%  Case 3: Serie RL with Sinusoidal source
%  For more information: https://electrica.ual.es/spacor/cases/index/3
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
R = 1.2;
L = 0.002;
I1 = 120;
w = 2*pi*50;
ii = I1*cos(w*t);
di = -w*I1*sin(w*t);
ddi = -w^2*I1*cos(w*t);
v = -R*ii-L*di;
dv = -R*di-L*ddi;
Rs = (di.*dv-v.*ddi)./(ii.*ddi-di.^2);
Ls = (v.*di-ii.*dv)./(ii.*ddi-di.^2);


subplot(2,1,1);
hold on;
grid on;
plot(t,v,'r','LineWidth',2);
plot(t,ii,'b','LineWidth',2);
text(0.005,140,'$v(t)$','interpreter','latex','fontsize',18);
text(0.018,170,'$i(t)$','interpreter','latex','fontsize',18);
set(gca,'FontSize',20,'TickLabelInterpreter','latex');

subplot(2,1,2);
hold on
grid on
plot(t,Rs,'r','LineWidth',2);
plot(t,500*Ls,'b','LineWidth',2);
ylim([0,1.5]);
xlabel('time (s)','interpreter','latex', 'FontWeight','bold','fontsize',16);
text(0.015,0.85,'$500L(t)$','interpreter','latex','fontsize',18);
text(0.010,1.3,'$R(t)$','interpreter','latex','fontsize',18);
set(gca,'FontSize',20,'TickLabelInterpreter','latex');

suptitle('Case 3: Serie RL with Sinusoidal source');
