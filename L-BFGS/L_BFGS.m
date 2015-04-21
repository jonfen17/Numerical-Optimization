%========================================================================
% Este programa lleva a cabo el m�todo
% L-BFGS
% para obtener la soluci�n al sistema de ecuaciones 
% A*x = b
% Curso An�lisis Aplicado
% Abril 2015
% Carlos Dioney Blanco Gonz�lez 
%========================================================================

function sol = L_BFGS( fun, x0, maxiter, m )

%========================================================================
% INPUT:
%  fun      - La funci�n a minimizar
%  x        - El punto inicial
% maxiter   - M�ximo de iteraciones
%  m        - Tama�o LBFGS
%
% OUTPUT:
%  sol     - Soluci�n del problema
%  i       - N�mero total de iteraciones
%========================================================================
tol = 1e-8;               % tolerancia
i   = 0;                  % iniciamos nuestro contador
f0  = feval(fun, x0);     % funci�n y gradiente en el punto inicial
g0  = diffgrad( x0, fun);
n   = length(x0);       
S   = zeros( n, m);       % matrices S, Y y Rho para calcular H*g
Y   = zeros( n, m);
Rho = zeros( m, 1);
fprintf(1,'   k        f            ||g||     alfa     nfg     curv  \n');
fprintf(1,'----------------------------------------------------------\n');
%========================================================================
g = g0;      x = x0;        p = -g;

[al, xn, f, gn, ~, nf]  = ...                       % BL con CFW
    linesch_sw( x, f0, g, p, fun, 0.0001, 0.9, 2); 
norma = norm(gn);
s     = xn - x;              % Calculamos S, Y y curvatura
y     = gn - g;
curv  = s' * y;
% Imprimimos y Actualzamos
fprintf(' %3i  %14.8e   %8.2e   %5.3f    %3i       %5.3f   \n', ...
                            i, f,  norma, al, nf, curv );
x = xn;
g = gn;
i = 1;

while (norma > tol) && (i < maxiter) && (curv > 0)
        
    h    = s'*y/ (y'*y);
    H    = h * ones( n, n);
    
    if   i > m
        
        % Actualizamos S, Y que entran a L-BFGS 
        % Notar que s�lo cambiamos una columna ( idx) en cada iteraci�n
        idx = mod(i,m);
        
        if idx == 0
            idx = m;
        end
        
        Rho(idx) = 1/curv;
        S(:,idx) = s;
        Y(:,idx) = y;
        p        = -DoubleLoop( g, H, S, Y, Rho, i, m);
        
    elseif i <= m
        
        % Actualizamos S,Y que entran a L-BFGS
        % En estos pasos se hace el BFGS normal
        Rho(i)  = 1/curv;
        S(:, i) = s;
        Y(:, i) = y;
        p       = -DoubleLoop( g, H, S(:,1:i), Y(:,1:i), Rho, i, m);
        
    end
    
    [al, xn, f, gn, ~, nf]  = ...                     % BL con CFW
        linesch_sw( x, f, g, p, fun, 0.0001, 0.9, 2);
    % Actualizamos 
    s     = xn - x;           % S y G para BFGS
    y     = gn - g;
    curv  = s' * y;    
    norma = norm(gn);
    x     = xn;
    g     = gn;
  
      fprintf(' %3i  %14.8e   %8.2e   %5.3f    %3i       %5.3f   \n', ...
                            i, f,  norma, al, nf, curv );
    i = i + 1;
end

sol = x;  % Regresamos la soluci�n