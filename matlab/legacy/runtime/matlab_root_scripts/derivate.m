function y_dot = derivate(y, dt,order)
    % DERIVATE calculates the numerical derivative of a time series.
    % 
    % y_dot = derivate(y, dt, N, order)
    % 
    % Inputs:
    %   y     - Input signal or time series.
    %   dt    - Time step or spacing between elements of y.
    %   N     - Number of points to use for differentiation, must be an odd integer >= 5.
    %   order - Order of the derivative to calculate.
    %
    % Output:
    %   y_dot - Calculated derivative of the input signal.
    %

    % Calculate coefficients
    b = computeDerivativeCoefficients(11, order);

    % Apply convolution
    y_dot_intermediate = conv(y, b, 'same');

    % Divide by dt^order
    y_dot = y_dot_intermediate / dt^order;
      
end


function b = computeDerivativeCoefficients(N, order)
    % Computes coefficients for robust numerical differentiation.
    %
    % b = computeDerivativeCoefficients(N, order)
    %
    % Inputs:
    %   N     - Odd integer >= 5.
    %   order - Positive integer specifying the derivative order.
    %
    % Output:
    %   b     - Coefficient array for the derivative filter.

    % Input checking
    assert(rem(N, 2) == 1 && N >= 5, 'N must be a positive odd integer >= 5.');
    assert(order > 0 && floor(order) == order, 'Order must be a positive integer.');

    % Calculate m and M
    m = (N - 3) / 2;
    M = (N - 1) / 2;
    k = M:-1:1;

    % Calculate coefficients for the given order
    % Construct the full filter coefficients
    c_k = (binomialCoefficient(2*m, m-k+1) - binomialCoefficient(2*m, m-k-1)) / 2^(2*m+1);
    switch order
        case 1           
            b = [c_k, 0, -c_k(end:-1:1)];
        case 2
            c_k = c_k.^2;
            b = 5/2*[c_k,0, -c_k(end:-1:1)-c_k,0, c_k(end:-1:1)];
        case 3
             b = 74/2500*[c_k,0,-c_k-c_k(end:-1:1)-c_k,0,c_k(end:-1:1)+c_k+c_k(end:-1:1),0,-c_k(end:-1:1)];       
        otherwise
            error('Higher order derivatives are not implemented for order %d', order);
    end

end

function coefficients = binomialCoefficient(n, k)
    % Computes binomial coefficients for each element in k using a scalar n.
    %
    % coefficients = binomialCoefficient(n, k)
    %
    % Inputs:
    %   n - A non-negative integer scalar.
    %   k - A vector of non-negative integers.
    %
    % Output:
    %   coefficients - A vector of binomial coefficients.

    % Preallocate coefficients
    coefficients = zeros(size(k));

    % Calculate coefficients for valid k values
    valid_idx = k >= 0 & k <= n;
    coefficients(valid_idx) = arrayfun(@(x) nchoosek(n, x), k(valid_idx));
end

function s = stirling_first_kind(n, k)
    % Calcula el coeficiente de Stirling de primera especie s(n, k)
    % Utiliza una tabla para evitar recalcular valores

    % Crear tabla (matriz) para almacenar valores
    s = zeros(n+1, k+1);
    
    % Inicializar la tabla
    for i = 0:n
        for j = 0:min(i, k)
            if j == 0 && i == 0
                s(i+1, j+1) = 1;
            elseif j == 0 || i == 0
                s(i+1, j+1) = 0;
            else
                s(i+1, j+1) = (i-1) * s(i, j+1) - s(i, j);
            end
        end
    end
    
    % El resultado es s(n, k)
    s = s(n+1, k+1);
end

function S = stirling_second_kind(n, k)
    % Calcula el coeficiente de Stirling de segunda especie S(n, k)
    % Utiliza una tabla para evitar recalcular valores

    % Crear tabla (matriz) para almacenar valores
    S = zeros(n+1, k+1);
    
    % Inicializar la tabla
    for i = 0:n
        for j = 0:min(i, k)
            if j == 0 && i == 0
                S(i+1, j+1) = 1;
            elseif j == 0
                S(i+1, j+1) = 0;
            else
                S(i+1, j+1) = j * S(i, j+1) + S(i, j);
            end
        end
    end
    
    % El resultado es S(n, k)
    S = S(n+1, k+1);
end


function coeff = generalizedBinomial(n, k)
    % Computes generalized binomial coefficient.
    % Works for integer and non-integer n and k values.
    %
    % coeff = generalizedBinomial(n, k)

    % Vectorize k if necessary
    if isscalar(k)
        k = [k];
    end

    % Preallocate result
    coeff = zeros(size(k));
    
    % Calculate coefficients using the Gamma function
    coeff(:) = gamma(n+1) ./ (gamma(k+1) .* gamma(n-k+1));
end