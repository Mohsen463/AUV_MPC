function [y,grady] = quadobj_LT(xLT,A2A_d,iii)
         y = (-1)^iii * A2A_d (1,:)*xLT ; 
  %% gradiant of the objective function (just take derivitave with respect to x)
if nargout > 1
     grady = A2A_d(1,:);
end