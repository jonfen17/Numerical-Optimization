% funci�n rosenbrock en R^2
% Esta funci�n alcanza su m�nimo en x = [ 1, 1 ]'

function [f, g, H] = rosen(x)
%function [f, g ] = rosen(x);
f = 100*( x(2) - x(1)^2 )^2  +  ( 1 - x(1) )^2 ;

g = [ -400*( x(2) - x(1)^2 )*x(1) - 2*( 1 - x(1) ) ;
       200*( x(2) - x(1)^2 )     ];

H = [ 1200*x(1)^2 - 400*x(2) + 2   -400*x(1) ;
       -400*x(1)                    200      ];