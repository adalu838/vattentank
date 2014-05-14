function [S,alarm] = decisioncalc(T,s0,s1)
% DECISIONCALC  Calculate the diagnosis S from thresholded residuals
%
%   Syntax: [S,alarm] = decisioncalc(T,s0,s1)
%
%   In:      T  The thresholded residuals. Each thresholded 
%               residual in separate columns.
%           s0  Binary table describing the decision when residuals is
%               below the threshold.               
%           s1  Binary table describing the decision when residuals is
%               above the threshold
%
%   Out: alarm  Signal indicating when a fault is detected
%            S  Binary matrix describing the diagnosis. A 1 on row i
%               column j indicates fault mode j is detected at time i.
%
%   Example:
%          Consider the decisionstructure
%
%          NF f1 f2 f3
%       T1 0  0  X  X
%       T2 0  X  0  1  
%       T3 0  X  X  0
%
%       s0 = [1 1 1 1;1 1 1 0;1 1 1 1]
%       s1 = [0 0 1 1;0 1 0 1;0 1 1 0]
%
%   Fault mode NF (No Fault) MUST be the first fault mode.

mT = size(T,2);
nT = size(T,1);

if (mT~=size(s0,1))|(mT~=size(s1,1))
  error('decstruc mismatch');
end
nF = size(s1,2);

si = zeros(nT,mT*nF);
for k=1:mT
  idx = find(T(:,k)==0);
  si(idx,(k-1)*nF+1:k*nF)=ones(length(idx),1)*s0(k,:);
  idx = find(T(:,k)==1);
  si(idx,(k-1)*nF+1:k*nF)=ones(length(idx),1)*s1(k,:);
end
S = ones(nT,nF);
for k=1:nF
  for l=1:mT
    S(:,k) = S(:,k).*si(:,(l-1)*nF+k);
  end;  
end
alarm = 1-S(:,1);
