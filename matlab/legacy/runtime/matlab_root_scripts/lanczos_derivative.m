% Script para procesar señales de tensión eléctrica con armónicos
% Utilizando filtro Savitzky-Golay para suavizado y derivación

% Parámetros del filtro Savitzky-Golay
orden_polinomio = 5;  % Orden del polinomio de ajuste
longitud_ventana = 25;  % Longitud de la ventana de filtrado (debe ser impar)

% Generar una señal de ejemplo con armónicos
Fs = 10000;  % Frecuencia de muestreo
t = 0:1/Fs:1;  % Vector de tiempo de 1 segundo

% Señal fundamental
fundamental = 50;  % Frecuencia fundamental (Hz)
amplitud_fundamental = 220;  % Amplitud de la señal fundamental (V)

% Armónicos
armonica_3 = 0.3 * amplitud_fundamental;  % Armónica de tercer orden
armonica_5 = 0.1 * amplitud_fundamental;  % Armónica de quinto orden

% Ruido
ruido = 5 * randn(size(t));  % Ruido gaussiano

% Construcción de la señal
x = amplitud_fundamental * sin(2*pi*fundamental*t) + ...
    armonica_3 * sin(2*pi*(3*fundamental)*t) + ...
    armonica_5 * sin(2*pi*(5*fundamental)*t) + ...
    ruido;

% Diseño del filtro Savitzky-Golay
[b, g] = sgolay(orden_polinomio, longitud_ventana);

% Calcular derivadas
derivadas = zeros(length(x), 4);

% Calcular señal suavizada y derivadas
for p = 0:3
    derivadas(:, p+1) = conv(x, factorial(p)/(-1/Fs)^p * g(:, p+1), 'same');
end

% Graficar resultados
figure;
subplot(2,2,1);
plot(t, x, 'DisplayName', 'Señal Original');
title('Señal Original');
xlabel('Tiempo (s)');
ylabel('Tensión (V)');
legend;

subplot(2,2,2);
plot(t, derivadas(:,1), 'DisplayName', 'Señal Suavizada');
title('Señal Suavizada');
xlabel('Tiempo (s)');
ylabel('Tensión (V)');
legend;

subplot(2,2,3);
plot(t, derivadas(:,2), 'DisplayName', 'Primera Derivada');
title('Primera Derivada');
xlabel('Tiempo (s)');
ylabel('dV/dt');
legend;

subplot(2,2,4);
plot(t, derivadas(:,3), 'DisplayName', 'Segunda Derivada');
title('Segunda Derivada');
xlabel('Tiempo (s)');
ylabel('d²V/dt²');
legend;

% Gráfico adicional para la tercera derivada
figure;
plot(t, derivadas(:,4), 'DisplayName', 'Tercera Derivada');
title('Tercera Derivada');
xlabel('Tiempo (s)');
ylabel('d³V/dt³');
legend;