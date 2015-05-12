%==========================================================================
% Este programa realiza el m�todo de Newton con 
% gradiente conjugado aplicado a la funci�n rosenbrock
% Lineal Avanzada       Mayo 2015
% Carlos Dioney Blanco Gonz�lez
% Stefano Molina Mart�nez
%==========================================================================
% Input:
% fun     - funci�n
% x       - punto inicial
% maxiter - m�ximo n�mero de iteraciones
%
% Output:
% x     - soluci�n
% i     - iteraciones
%--------------------------------------------------------------------------

function [x] = NewtonGC(fun, x, maxiter)
% Definimos par�metros 
n     = length(x);
c1    = 10e-4;
c2    = .9;
tol   = 1.0e-6;
tol2  = 1.0e-12;
maxCG = 2*n;

[f1, g, H] = feval(fun, x);
norma      = norm(g);

i = 0;
fprintf(1,'    iter      alpha        ||grad(f)||         f        \n\n');

while (norma > tol) && (i < maxiter)
    
    % Obtenemos direcci�n por GC
    p = gradconj( H, -g, tol, maxCG, tol2);
    
    % Verificamos tama�o del paso con Condiciones de Wolfe
    alpha = 1;
    [f2, g2, ~] = feval( fun, (x + alpha*p));
    
    while  ( f2 > f1 + c1*alpha*g'*p ) || ( abs(p'*g2) > c2*abs(p'*g) )
        % Bisecci�n
        alpha = alpha/2;
        [f2, g2, ~] = feval(fun, (x + alpha*p));
        
    end
    
    % Actualizamos a nuevo punto, gradiente y Hessiana.
     x         = x + alpha*p;
    [f1, g, H] = feval(fun, x);
    
    i     = i + 1;
    norma = norm(g);
    
    fprintf(1,'    %3.0f     %5f     %5.5f        %5.5f    \n',i,alpha,norma,f1); 
end
    if (i == maxiter)
    fprintf(1, 'El proceso alcanz� el m�ximo n�mero de iteraciones. \n\n' );
    else
    fprintf(1, 'Se realiz� el proceso con �xito. \n\n');
    end
disp('Se alcanz� el m�nimo en:');
x
end