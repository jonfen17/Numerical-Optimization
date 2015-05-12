%==========================================================================
% Este programa realiza el método de Newton con 
% gradiente conjugado aplicado a la función rosenbrock
% Lineal Avanzada       Mayo 2015
% Carlos Dioney Blanco González
% Stefano Molina Martínez
%==========================================================================
% Input:
% fun     - función
% x       - punto inicial
% maxiter - máximo número de iteraciones
%
% Output:
% x     - solución
% i     - iteraciones
%--------------------------------------------------------------------------

function [x] = NewtonGC(fun, x, maxiter)
% Definimos parámetros 
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
    
    % Obtenemos dirección por GC
    p = gradconj( H, -g, tol, maxCG, tol2);
    
    % Verificamos tamaño del paso con Condiciones de Wolfe
    alpha = 1;
    [f2, g2, ~] = feval( fun, (x + alpha*p));
    
    while  ( f2 > f1 + c1*alpha*g'*p ) || ( abs(p'*g2) > c2*abs(p'*g) )
        % Bisección
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
    fprintf(1, 'El proceso alcanzó el máximo número de iteraciones. \n\n' );
    else
    fprintf(1, 'Se realizó el proceso con éxito. \n\n');
    end
disp('Se alcanzó el mínimo en:');
x
end