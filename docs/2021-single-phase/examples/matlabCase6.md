# Matlab simulation for case 6
## Parallel RLC with switching capacitor and resistor
<img src="https://electrica.ual.es/spacor/images/spacorcasesix_simulation.png" width="100%">

Matlab Code
Herein, you can find the code for simulation purposes.

```matlab                   
%**************************************************************************
%  SPACOR Space Vectors for Parameter Identification of Non-Linear Circuits
%  https://electrica.ual.es/spacor
%
%  Case 6: Parallel RLC with switching capacitor and resistor
%  For more information: https://electrica.ual.es/spacor/cases/index/6
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
Vh = 12;
h = 7;
G1 = 1/50;
G2 = 1/77;
G_aux = ones(1,length(t));
G = [G1*G_aux(1:end/2) G2*G_aux(end/2+1:end)];
C1 = 1e-5;
C2 = 3e-6;
C_aux = ones(1,length(t));
C = [C1*C_aux(1:350) C2*C_aux(351:772) C1*C_aux(773:1058) (C1+C2)*C_aux(1059:end)];
L1 = 0.7;
 
v = V1*sin(w*t)+Vh*sin(h*w*t);
dv = w*V1*cos(w*t)+(w*h)*Vh*cos(h*w*t);
ddv = -w^2*V1*sin(w*t)-(w*h)^2*Vh*sin(h*w*t);
dddv = -w^3*V1*cos(w*t)-(w*h)^3*Vh*cos(h*w*t);
ddddv = w^4*V1*sin(w*t)+(w*h)^4*Vh*sin(h*w*t);

intv = -V1/w*cos(w*t)-Vh/(h*w)*cos(h*w*t);
iL = 1/L1*intv;
diL = 1/L1*v;
ddiL = 1/L1*dv;
dddiL = 1/L1*ddv;

iR = G.*v;
diR = G.*dv;
ddiR = G.*ddv;
dddiR = G.*dddv;

iC = C.*dv;
diC = C.*ddv;
ddiC = C.*dddv;
dddiC = C.*ddddv;

ii = iR+iL+iC;
di = diR+diL+diC;
ddi = ddiR+ddiL+ddiC;
dddi = dddiR+dddiL+dddiC;

Gp = (dddi.*(v.*dddv-dv.*ddv)-ddddv.*(v.*ddi-dv.*di)+ddv.*(ddv.*ddi-dddv.*di))...
     ./(ddddv.*(dv.^2-ddv.*v)-ddv.*(dv.*dddv-ddv.^2)+dddv.*(v.*dddv-dv.*ddv));
Gamma = -(dddi.*(dv.*dddv-ddv.^2)-ddddv.*(dv.*ddi-ddv.*di)+dddv.*(ddv.*ddi-dddv.*di))...
        ./(ddddv.*(dv.^2-ddv.*v)-ddv.*(dv.*dddv-ddv.^2)+dddv.*(v.*dddv-dv.*ddv));
Cap = (dddi.*(dv.^2-ddv.*v)-ddv.*(dv.*ddi-ddv.*di)+dddv.*(v.*ddi-dv.*di))...
      ./(ddddv.*(dv.^2-ddv.*v)-ddv.*(dv.*dddv-ddv.^2)+dddv.*(v.*dddv-dv.*ddv));
Gp = filloutliers(Gp,'nearest','movmedian',100);
Gamma = filloutliers(Gamma,'nearest','movmedian',100);
Cap = filloutliers(Cap,'nearest','movmedian',100);

subplot(3,1,1);
hold on;
grid on;
grid(gca,'minor');
plot(t,iC,'b','LineWidth',2);
legend({'$i_C(t)$'},'Location','southwest','interpreter','latex');
set(gca,'FontSize',16,'TickLabelInterpreter','latex');

subplot(3,1,2)
hold on;
grid on;
grid(gca,'minor');
plot(t,v,'r','LineWidth',2);
plot(t,40*ii,'b','LineWidth',2);
legend({'$v(t)$','$40i(t)$'},'Location','southwest','interpreter','latex');
set(gca,'FontSize',16,'TickLabelInterpreter','latex');

subplot(3,1,3)
hold on;
grid on;
grid(gca,'minor');
plot(t,100*Gp,'r','LineWidth',2);
plot(t,Gamma,'b','LineWidth',2);
plot(t,0.5e5*Cap,'c','LineWidth',2);
ylim([0,2.2]);
xlabel('time (s)','interpreter','latex', 'FontWeight','bold','fontsize',16);
legend({'$100G(t)$','$\Gamma(t)$','$0.5\times10^5C(t)$'},'Location','northeast','interpreter','latex');
set(gca,'FontSize',16,'TickLabelInterpreter','latex');

suptitle('Case 6: Parallel RLC with switching capacitor and resistor');
