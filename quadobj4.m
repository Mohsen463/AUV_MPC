function [y,grady] = quadobj4(x4,ii,A2A_d)
         y = - A2A_d (1,:,ii)*x4 ; 
  %% gradiant of the objective function (just take derivitave with respect to x)
if nargout > 1
     grady = A2A_d(1,:,ii);
end