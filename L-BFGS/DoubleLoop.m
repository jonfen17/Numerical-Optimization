%-------------------------------------------------------------------------
% Carlos Dioney Blanco González     131370
% Curso: Análsis Aplicado           Primavera 2015
% L-BFGS recursión de doble ciclo
%
% Este programa realiza la multiplicación de una matriz H y un vector g
% correspondinte a la aproximación de la hessiana por BFGS y el gradiente
% respectivamente. 
%
% Input:
% g - Gradiente en la k-ésima iteración
% H - Aproximación inicial de la hessiana
% S - Matriz con columnas dadas por  S(:,i) = X_i+1 - X_i   i = k-m ... k-1
% Y - Matriz con columnas dadas por  Y(:,i) = g_i+1 - g_i   i = k-m ... k-1 
% Output:
% r - El producto de H*g
%--------------------------------------------------------------------------

function r = DoubleLoop( g, H, S, Y)

    [ n, m] = size( S );
    cero = zeros( m, 1);  
    rho	    = cero;  	alpha   = cero; 	beta    = cero;

    % Calculamos nuestras rho en un vector de k*1
    for i = 1:m
       rho(i,1) = 1/(Y(:,i)'*S(:,i)); 
    end 

    % Realizamos el primer loop
    for i = m:-1:1

        alpha(i) = rho(i)* S(:,i)'*g;
        g        = g - alpha(i)* Y(:,i);

    end

    r = H*g;
    
    % Realizamos el segundo loop
    for i = 1:m

        beta(i) = rho(i)* Y(:,i)'*r;
        r       = r + S(:,i)*(alpha(i)- beta(i));
        
    end

end