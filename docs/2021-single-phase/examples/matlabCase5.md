# Matlab simulation for case 5
## Parallel RL with piecewise non-linear inductor
<img src="https://electrica.ual.es/spacor/images/spacorcasefive_simulation.png" width="100%">

Matlab Code
Herein, you can find the code for simulation purposes.

```matlab                   
%**************************************************************************
%  SPACOR Space Vectors for Parameter Identification of Non-Linear Circuits
%  https://electrica.ual.es/spacor
%
%  Case 5: Parallel RL with piecewise non-linear inductor
%  For more information: https://electrica.ual.es/spacor/cases/index/5
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
DT = 0.00001;
t = [0:DT:0.02-DT];
V1 = sqrt(2)*120;
w = 2*pi*50;
G1 = 1/0.5;
G2 = 1;
G_aux = ones(1,length(t));
G = [G1*G_aux(1:end/2) G2*G_aux(end/2+1:end)];
v = V1*cos(w*t);
dv = -w*V1*sin(w*t);
ddv = -w^2*V1*cos(w*t);
intv = V1/w*sin(w*t);
fs = 0.4;
L1 = 0.007;
L2 = 0.001;

idx1=find(intvfs-0.0001);
idx2=find(intv>-fs-0.0009&intv<-fs+0.0001);
sat=sort([idx1 idx2]);

iL1 = [1/L1*intv(1:sat(1)) 1/L2*intv(sat(1)+1:sat(2))+0.4*(1/L1-1/L2)...
    1/L1*intv(sat(2)+1:(sat(3))) 1/L2*intv((sat(3)+1):sat(4))-0.4*(1/L1-1/L2)...
    1/L1*intv(sat(4)+1:length(t))];
diL1 = [1/L1*v(1:sat(1)) 1/L2*v(sat(1)+1:sat(2))...
    1/L1*v(sat(2)+1:(sat(3))) 1/L2*v((sat(3)+1):sat(4))...
    1/L1*v(sat(4)+1:length(t))];
ddiL1 = [1/L1*dv(1:sat(1)) 1/L2*dv(sat(1)+1:sat(2))...
    1/L1*dv(sat(2)+1:(sat(3))) 1/L2*dv((sat(3)+1):sat(4))...
    1/L1*dv(sat(4)+1:length(t))];

iR = G.*v;
diR = G.*dv;
ddiR = G.*ddv;

ii = iR+iL1;
di = diR+diL1;
ddi = ddiR+ddiL1;
GG = (v.*ddi-di.*dv)./(v.*ddv-dv.^2);
GaGa = -(dv.*ddi-di.*ddv)./(v.*ddv-dv.^2);

subplot(3,1,1);
hold on;
grid on;
grid(gca,'minor');
plot(t,iL1,'b','LineWidth',2);
legend({'$i_L(t)$'},'Location','southwest','interpreter','latex');
set(gca,'FontSize',16,'TickLabelInterpreter','latex');

subplot(3,1,2);
hold on;
grid on;
grid(gca,'minor');
plot(t,2*v,'r','LineWidth',2);
plot(t,ii,'b','LineWidth',2);
legend({'$2v(t)$','$i(t)$'},'Location','southwest','interpreter','latex');
set(gca,'FontSize',16,'TickLabelInterpreter','latex');

subplot(3,1,3);
hold on;
grid on;
grid(gca,'minor');
plot(t,500*GG,'r','LineWidth',2);
plot(t,GaGa,'b','LineWidth',2);
ylim([0,1100]);
xlabel('time (s)','interpreter','latex', 'FontWeight','bold','fontsize',16);
legend({'$500G(t)$','$\Gamma(t)$'},'Location','northwest','interpreter','latex');
set(gca,'FontSize',16,'TickLabelInterpreter','latex');

suptitle('Case 5: Parallel RL with piecewise non-linear inductor');

