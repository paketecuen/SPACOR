# Matlab simulation for case 7
## Series RLC with switching capacitor and resistor
<img src="https://electrica.ual.es/spacor/images/spacorcaseseven_simulation.png" width="100%">

Matlab Code
Herein, you can find the code for simulation purposes.

```matlab                   
%**************************************************************************
%  SPACOR Space Vectors for Parameter Identification of Non-Linear Circuits
%  https://electrica.ual.es/spacor
%
%  Case 7: Series RLC with switching capacitor and resistor
%  For more information: https://electrica.ual.es/spacor/cases/index/7
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
I1 = sqrt(2)*10;
w = 2*pi*50;
Ih = 2;
h = 7;
R1 = 5;
R2 = 15;
R_aux = ones(1,length(t));
R = [R1*R_aux(1:end/2) R2*R_aux(end/2+1:end)];
C1 = 1e-1;
C2 = 3e-2;
C_aux = ones(1,length(t));
C = [C1*C_aux(1:350) C2*C_aux(351:772) C1*C_aux(773:1058) (C1+C2)*C_aux(1059:end)];
L1 = 0.007;

ii = I1*sin(w*t)+Ih*sin(h*w*t);
di = w*I1*cos(w*t)+(w*h)*Ih*cos(h*w*t);
ddi = -w^2*I1*sin(w*t)-(w*h)^2*Ih*sin(h*w*t);
dddi = -w^3*I1*cos(w*t)-(w*h)^3*Ih*cos(h*w*t);
ddddi = w^4*I1*sin(w*t)+(w*h)^4*Ih*sin(h*w*t);
inti = -I1/w*cos(w*t)-Ih/(h*w)*cos(h*w*t);

vL = L1*di;
dvL = L1*ddi;
ddvL = L1*dddi;
dddvL = L1*ddddi;

vR = R.*ii;
dvR = R.*di;
ddvR = R.*ddi;
dddvR = R.*dddi;

vC = 1./C.*inti;
dvC = 1./C.*ii;
ddvC = 1./C.*di;
dddvC = 1./C.*ddi;

v = vR+vL+vC;
dv = dvR+dvL+dvC;
ddv = ddvR+ddvL+ddvC;
dddv = dddvR+dddvL+dddvC;


Rs = (dddv.*(ii.*dddi-di.*ddi)-ddddi.*(ii.*ddv-di.*dv)+ddi.*(ddi.*ddv-dddi.*dv))...
     ./(ddddi.*(di.^2-ddi.*ii)-ddi.*(di.*dddi-ddi.^2)+dddi.*(ii.*dddi-di.*ddi));
S = -(dddv.*(di.*dddi-ddi.^2)-ddddi.*(di.*ddv-ddi.*dv)+dddi.*(ddi.*ddv-dddi.*dv))...
    ./(ddddi.*(di.^2-ddi.*ii)-ddi.*(di.*dddi-ddi.^2)+dddi.*(ii.*dddi-di.*ddi));
Ls = (dddv.*(di.^2-ddi.*ii)-ddi.*(di.*ddv-ddi.*dv)+dddi.*(ii.*ddv-di.*dv))...
     ./(ddddi.*(di.^2-ddi.*ii)-ddi.*(di.*dddi-ddi.^2)+dddi.*(ii.*dddi-di.*ddi));
Rs = filloutliers(Rs,'nearest','movmedian',100);
S = filloutliers(S,'nearest','movmedian',100);
Ls = filloutliers(Ls,'nearest','movmedian',100); 
Rs = filloutliers(Rs,'nearest','movmedian',100);
S = filloutliers(S,'nearest','movmedian',100);
Ls = filloutliers(Ls,'nearest','movmedian',100); 

subplot(3,1,1);
hold on;
grid on;
grid(gca,'minor');
plot(t,100*vC,'b','LineWidth',2);
plot(t,4*vL,'r','LineWidth',2);
plot(t,vR,'c','LineWidth',2);
legend({'$100v_C(t)$','$4v_L(t)$','$v_R(t)$'},'Location','southwest','interpreter','latex');
set(gca,'FontSize',16,'TickLabelInterpreter','latex');

subplot(3,1,2);
hold on;
grid on;
grid(gca,'minor');
plot(t,v,'r','LineWidth',2);
plot(t,10*ii,'b','LineWidth',2);
legend({'$v(t)$','$10i(t)$'},'Location','southwest','interpreter','latex');
set(gca,'FontSize',16,'TickLabelInterpreter','latex');

subplot(3,1,3);
hold on;
grid on;
grid(gca,'minor');
plot(t,Rs,'r','LineWidth',2);
plot(t,4e3*Ls,'b','LineWidth',2);
plot(t,S,'c','LineWidth',2);
xlabel('time (s)','interpreter','latex', 'FontWeight','bold','fontsize',16) ;                      
legend({'$R(t)$','$4\times10^3L(t)$','$S(t)$'},'Location','northeast','interpreter','latex');
set(gca,'FontSize',16,'TickLabelInterpreter','latex');

suptitle('Case 7: Series RLC with switching capacitor and resistor');
