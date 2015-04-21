%-------------------------------------------------------------------------
% Carlos Dioney Blanco Gonz�lez     131370
% Curso: An�lsis Aplicado           Primavera 2015
% L-BFGS recursi�n de doble ciclo
%
% Este programa realiza la multiplicaci�n de una matriz H y un vector g
% correspondinte a la aproximaci�n de la hessiana por BFGS y el gradiente
% respectivamente. 
%
% Input:
% g   - Gradiente en la k-�sima iteraci�n
% H   - Aproximaci�n inicial de la hessiana
% S   - Matriz c/ columnas  S(:,i) = X_i+1 - X_i   i = k-m ... k-1
% Y   - Matriz c/ columnas  Y(:,i) = g_i+1 - g_i   i = k-m ... k-1
% rho - Vector con entradas 1 / (y' * s); de las m �ltimas iteraciones
% k   - Iteraci�n en BFGS
% m   - Columnas matriz S, Y
% Output:
% r - El producto de H*g
%--------------------------------------------------------------------------

function r = DoubleLoop( g, H, S, Y, rho, k, m)
    
    m       = min( k, m);
    cero    = zeros( 1, m);  
    alpha   = cero;
    idx     = backIndex( k, m); % Apuntador hacia atr�s

    % Realizamos el primer loop
    for i = 1:m
        j        = idx(i);
        alpha(j) = rho(j)* S(:,j)'*g;
        g        = g - alpha(j)* Y(:,j);

    end

    r   = H*g;
    idx = forwIndex( k, m); % Apuntador hacia adelante
    
    % Realizamos el segundo loop
    for i = 1:m
        j    = idx(i);
        beta = rho(j)* Y(:,j)'*r;
        r    = r + S(:,j)*(alpha(j)- beta);
        
    end

end