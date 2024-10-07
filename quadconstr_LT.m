function [y,yeq,grady,gradyeq] = quadconstr_LT(xLT,Nodes_num,a2a_num)
  N = Nodes_num;
  Num_ineq_const = 2*a2a_num + 2*(Nodes_num);
  y = zeros(1, Num_ineq_const);
  grady = zeros(a2a_num,Num_ineq_const);
  
  for jj = 1 : a2a_num
      cons_aux1 = zeros(1,a2a_num); cons_aux1(1,jj) = 1;
          y(1,jj) = cons_aux1 * xLT - 1;
      grady(:,jj) = cons_aux1';
  end

  for jj = 1 : a2a_num
      cons_aux2 = zeros(1,a2a_num); cons_aux2(1,jj) = -1;
          y(  a2a_num + jj) = cons_aux2 * xLT ;
      grady(:,a2a_num + jj) = cons_aux2';
  end 

  kk = 1;
  for jj = 1 : N  
      cons_aux3 = zeros(1,a2a_num);
      cons_aux3(1,kk:kk+(N-jj-1)) = 1;
      if jj > 1
          for ii = 1:jj-1
              if ii == 1
                 zz = jj-1 ;
              else
                  zz = zz + (N-ii)  ;
              end
              cons_aux3(1,zz) = 1;
          end
      end
           y(  2*a2a_num + (jj-1)*2+1) =  cons_aux3 * xLT - 2;
           y(  2*a2a_num + (jj-1)*2+2) = -cons_aux3 * xLT + 1;
       grady(:,2*a2a_num + (jj-1)*2+1) =  cons_aux3';
       grady(:,2*a2a_num + (jj-1)*2+2) = -cons_aux3';
       kk = kk+(N-jj);
  end 
 
cons_aux4 = ones(1,a2a_num);
yeq = cons_aux4 * xLT - (Nodes_num-1);
gradyeq = cons_aux4;
