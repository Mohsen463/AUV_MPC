function [y,grady] = quadobj3(x3,ii,USV_AUV_d)
         y =  USV_AUV_d (ii,:)*x3 ; 
  %% gradiant of the objective function (just take derivitave with respect to x)
if nargout > 1
     grady = USV_AUV_d (ii,:);
end