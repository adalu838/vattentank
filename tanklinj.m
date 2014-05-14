function Gsys=tanklinj(h10, d)
% Anrop: Gsys=tanklinj(h10,d)
%
%
% Skapar linjärisering Gsys av vattentankssystemet för vattennivå h1=h10.
% d är en vektor med modellparametrar
%

% parametrar enligt (7) i labb-PM.
d1 = d(1); d2 = d(2); d3 = d(3); d4 = d(4); d5 = d(5); d6 = d(6); 

% Bestämmer jämviktspunkt i tank 2
h20 = (d3)^2/(d4)^2*h10;


% dot{x} = f(x) = f(x0) + J(x0)(x-x0) + ... 
% y(x) = g(x0) + c*(x-x0) + ...
% där J(x0) = a
a=[-d2/(2*sqrt(h10)) 0; ...
    d3/(2*sqrt(h10)) -d4/(2*sqrt(h20))];

b= [d1; 0];

c=[1 0;... 
   0 1;...
   d5/(2*sqrt(h10)) 0;...
   0 d6/(2*sqrt(h20))];

d=[0;0;0;0];

Gsys=ss(a,b,c,d);