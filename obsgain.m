function K=obsgain(A,C,p)
% OBSGAIN  Pole placement technique to compute observer gain
%
%   K=obsgain(A,C,P) computes observer gain K such that 
%   the eigenvalues of matrix A-KC match the specified 
%   locations P. 

%
%   Method limitation: no eigenvalue can have a multiplicity
%   greater than the number of measurements.
%
%
%   See also: place, acker, lqe

n = size(A,1);
disp('Checking observability...');
if(rank(obsv(A,C))<n)
  error('System not observable');  
end;
disp('Observability ok.');

ps = unique(p);
my = 1;
for i=1:length(ps)
  my = max(my,sum(p==ps(i)));
end
if my>rank(C)
  error('Can''t place poles with multiplicity greater than rank(C).')
end;

K = place(A',C',p)';

