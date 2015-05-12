%==========================================================================
% Este programa lleva a cabo el m�todo del GC computacional
% Lineal Avanzada       Mayo 2015
% Carlos Dioney Blanco Gonz�lez
% Stefano Molina Mart�nez
%==========================================================================
% Input:
% H       - matriz
% g       - Cosntantes
% tol     - tolerancia del m�todo
% tol2    - para evitar curvatura negativa
%
% Output:
% x     - soluci�n
% i     - iteraciones
%--------------------------------------------------------------------------
% Ejemplo 
% A = [ 1 2 3 4 ; 2 4 6 7 ; 3 6 1 3 ; 4 7 3 8],   b = [ 1 3 5 6]'
%--------------------------------------------------------------------------
function [x,y] = gradconj(H, g, tol)
tic
% Condiciones iniciales para el gradiente conjugado
maxiter = 1000;
n = length(g);
x = zeros(n,1);
r = g + H*x;
d = r;
k = 0;

fprintf(1,'    iter         ||grad(r)||         \n\n');
% Implementamos GC para obtener una nueva direcci�n para Newton
while ( norm(r) > tol && k < maxiter)
    
    aux  = H * d; 
    daux = d'*aux;
    
    % GC computacional
    alpha  = ( r'*r )/( daux );
    x      = x + alpha*d;
    rvieja = r;
    r      = r - alpha*aux;
    beta   = (r'*r) / ( rvieja'*rvieja );
    dvieja = d;
    d      = r + beta*dvieja;
    k      = k + 1;
    
    % Imprimimos norma del residuo e iteraci�n
    fprintf(1,'    %3.0f         %5.5f       \n', k, norm(r) ); 

end
    % Imprimimos resultados
    if (k == maxiter)
        fprintf(1, 'El proceso alcanz� el m�ximo n�mero de iteraciones. \n\n' );
    else
        fprintf(1, 'Se realiz� el proceso con �xito. \n\n');
    end
    
disp('Soluci�n al problema:');
x
y = toc;
end