# Matlab simulation for case 1
## Parallel RL load with Sinusoidal source
<img src="https://electrica.ual.es/spacor/images/spacorcaseone_simulation.png" width="100%">

Matlab Code
Herein, you can find the code for simulation purposes.

```matlab                   
%**************************************************************************
%  SPACOR Space Vectors for Parameter Identification of Non-Linear Circuits
%  https://electrica.ual.es/spacor
%
%  Case 1: Parallel RL load with Sinusoidal source
%  For more information: https://electrica.ual.es/spacor/cases/index/1
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
DT = 0.00000001;
t = 0:DT:0.02-DT;
G1 = 1/1;
G2 = 1/2;
G_aux = ones(1,length(t));
G = [G1*G_aux(1:1000000) G2*G_aux(1000001:end)];
Gam1 = 1/0.0025;
Gam2 = 1/0.007;
Gam_aux = ones(1,length(t));
Gam = [Gam1*Gam_aux(1:1000000) Gam2*Gam_aux(1000001:end)];
V1 = 120;
w = 2*pi*50;
v = V1*sin(w*t);
dv = w*V1*cos(w*t);
ddv = -w^2*V1*sin(w*t);
dddv = -w^3*V1*cos(w*t);
iv = -V1/w*cos(w*t);
ii = G.*v+Gam.*iv;
di = G.*dv+Gam.*v;
dG = ((ii.*v)-(iv.*di))./(v.^2-iv.*dv);
dGam = ((v.*di)-(ii.*dv))./(v.^2-iv.*dv);

subplot(2,1,1);
hold on;
grid on;
plot(t,v,'r','LineWidth',2);
plot(t,ii,'b','LineWidth',2);
text(0.002,140,'$v(t)$','interpreter','latex','fontsize',18);
text(0.010,230,'$i(t)$','interpreter','latex','fontsize',18);
set(gca,'FontSize',20,'TickLabelInterpreter','latex');

subplot(2,1,2);
hold on;
grid on;
plot(t,500*dG,'r','LineWidth',2);
plot(t,dGam,'b','LineWidth',2);
xlabel('time (s)','interpreter','latex', 'FontWeight','bold','fontsize',12);
text(0.016,300,'$500G(t)$','interpreter','latex','fontsize',18);
text(0.004,420,'$\Gamma(t)$','interpreter','latex','fontsize',18);
set(gca,'FontSize',20,'TickLabelInterpreter','latex');

suptitle('Case 1: Parallel RL load with Sinusoidal source');
