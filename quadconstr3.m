function [y,yeq,grady,gradyeq] = quadconstr3(x3,auv_n)
  Num_ineq_const = 2*auv_n;
  y = zeros(1, Num_ineq_const);
  grady = zeros(auv_n,Num_ineq_const);
  
  for jj = 1 : auv_n
      cons_aux1 = zeros(1,auv_n); cons_aux1(1,jj) = 1;
          y(1,jj) = cons_aux1 * x3 - 1;
      grady(:,jj) = cons_aux1';
  end

  for jj = 1 : auv_n
      cons_aux2 = zeros(1,auv_n); cons_aux2(1,jj) = -1;
          y(  auv_n + jj) = cons_aux2 * x3 ;
      grady(:,auv_n + jj) = cons_aux2';
  end  
 
cons_aux3 = ones(1,auv_n);
yeq = cons_aux3 * x3 - 1;
gradyeq = cons_aux3;
