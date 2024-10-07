function [y,yeq,grady,gradyeq] = quadconstr5(x5,delta,ds,A2V_tol,Vmax,P_0,Z_floor,usv_head,auv_n,s_no,NumSlot,s_no_mpc,A_g_MPC_all,B_g_MPC_all,Z_AUV_MPC_all,P_usv_MPC,sig_mat,sig_ij)
% modeling the inequality constraints  like y( 1 , 1:Num_ineq_const)
% modeling the grad of inequality constraints  like y( s_no_mpc , 1:Num_ineq_const) and for each column(constraints) just take derivitave of each constraint in y with respect to x
    C1_4 = 1;         C1_6 = 1;              C1_7 = 1;      C1_14 = 1;      C1_12 = 1;    C5_2 = 1;      C5_3 = 1;% if constrant X (i.e., C1_X) is taken into acount it is 1, otherwise 0.
 % negative z;   dis(auv(K)-Usv(K))<ds;        speed         Heading          floor        A2U graph      A2A graph

Num_C1_4 =  C1_4 * auv_n*NumSlot;       Num_C1_6 = C1_6 * auv_n;             Num_C1_7 = C1_7 * (auv_n+1)*NumSlot;  % number of inequality constraints each C1_X imposes
Num_C1_14 = C1_14 * 4*NumSlot;            Num_C1_12 = C1_12 * auv_n*NumSlot;     Num_C5_2 = C5_2 * 2 * NumSlot;   Num_C5_3 = C5_3 * 2*NumSlot*(auv_n-1); 
Num_ineq_const = Num_C1_4  +  Num_C1_6  +  Num_C1_7  +  Num_C1_14  +  Num_C1_12 + Num_C5_2 + Num_C5_3;   % number of all inequality constraints  

y = zeros(1, Num_ineq_const);
grady = zeros(s_no_mpc,Num_ineq_const);

A_P0 = A_g_MPC_all*P_0;  B_delt = delta*B_g_MPC_all;
%% Constraint << C1_4 >> and its grad : inequality constraint associated with z should be negative (underwater sea level) 
 % number of inequality constraints Num_C1_4: auv_n * NumSlot
 if C1_4 == 1
        y(  1: auv_n * NumSlot) = (Z_AUV_MPC_all * ( A_P0 + B_delt*x5 ))'; 
    grady(:,1: auv_n * NumSlot) = (Z_AUV_MPC_all * B_delt )';
 end   

 %% Constraint C1_6 and its grad : unequality constraint associated with distances of AUVs from USV at the final stage be lower than a given distance ds based on the sonar modem/sensor rnge
 % number of inequality constraints Num_C1_6 = auv_n 
 if C1_6 == 1
    P_uav_individual = zeros(3,s_no_mpc);   % C1_6
    for jj = 1 : auv_n
        P_uav_individual(1:3 , (NumSlot-1)*s_no+(jj-1)*3+4 : (NumSlot-1)*s_no+(jj-1)*3+6) = eye(3);
            y(  Num_C1_4 + jj) = 1*  ((P_usv_MPC - P_uav_individual) * ( A_P0 + B_delt*x5))' * ((P_usv_MPC - P_uav_individual) * ( A_P0 + B_delt*x5)) - ds^2 ;
        grady(:,Num_C1_4 + jj) = 2*( ((P_usv_MPC - P_uav_individual)*B_delt)'  *  (  ((P_usv_MPC - P_uav_individual)*A_P0) + ((P_usv_MPC - P_uav_individual)*(B_delt*x5))  )  );
        P_uav_individual = zeros(3,s_no_mpc);
    end  
 end 
 %% Constraint C1_7 and its grad : unequality constraint associated with maximum input effort (speed limit)
  % number of inequality constraints Num_C1_7 = (auv_n+1)*NumSlot
if C1_7 == 1
   V_individual = zeros(3,s_no_mpc);   % C1_7 
  for ii = 1 : NumSlot
      for jj = 1 : auv_n + 1
          V_individual(1:3 , (ii-1)*s_no+(jj-1)*3+1 : (ii-1)*s_no+(jj-1)*3+3) = eye(3);
              y(  Num_C1_4 + Num_C1_6 + (ii-1)*(auv_n + 1)+jj ) =  (V_individual * x5)' * (V_individual * x5) - Vmax^2  ;       %%% NumSlot must be (auv_n + 1)
          grady(:,Num_C1_4 + Num_C1_6 + (ii-1)*(auv_n + 1)+jj ) =  (V_individual)'     * (V_individual * x5) ;  % gradiant     %%% NumSlot must be (auv_n + 1)
          V_individual = zeros(3,s_no_mpc); 
      end 
  end
end 

%% Constraint C1_14 and its grad : unequality constraint associated with heading of USV
  % number of inequality constraints Num_C1_14 = 4*NumSlot
if C1_14 == 1
   V_usv_x = zeros(1,s_no_mpc);V_usv_y = zeros(1,s_no_mpc);   % C1_14 
  for ii = 1 : NumSlot
      if ii == 1
         V_usv_x(1 , (ii-1)*s_no+1) = 1;  V_usv_y(1 , (ii-1)*s_no+2) = 1;
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+1 ) = (V_usv_x * x5)  -  (Vmax*cos(usv_head)) ; % for the first step time specify max of input speed in x & y direction based on the initial heading 
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+2 ) = (V_usv_y * x5)  -  (Vmax*sin(usv_head)) ; % for the first step time specify min of input speed in x & y direction be zero (to be positive) to keep in the first quarter of cartesian coordinates 
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+3 ) = -(V_usv_x * x5) +  (Vmax*cos(1.5*usv_head)) ; 
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+4 ) = -(V_usv_y * x5) +  (Vmax*sin(usv_head/2)) ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+1) =  (V_usv_x )'  ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+2) =  (V_usv_y )'  ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+3) = -(V_usv_x )'  ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+4) = -(V_usv_y )'  ;
         V_usv_x = zeros(1,s_no_mpc); V_usv_y = zeros(1,s_no_mpc); 
      else
         V_usv_x(1 , (ii-1)*s_no+1) = 1;  V_usv_x(1 , (ii-2)*s_no+1) = -(1+(1/Vmax)); % for modeling Vx[k] - Vx[k-1] < Vx[k-1]/Vmax
         V_usv_y(1 , (ii-1)*s_no+2) = 1;  V_usv_y(1 , (ii-2)*s_no+2) = -(1+(1/Vmax));
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+1) =  (V_usv_x * x5) ;
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+2) =  (V_usv_y * x5) ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+1) =  (V_usv_x )' ;        % gradiant
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+2) =  (V_usv_y )' ;
         V_usv_x = zeros(1,s_no_mpc); V_usv_y = zeros(1,s_no_mpc);

         V_usv_x(1 , (ii-1)*s_no+1) = -1;  V_usv_x(1 , (ii-2)*s_no+1) = (1-(1/Vmax));  % for modeling Vx[k] - Vx[k-1] > -Vx[k-1]/Vmax
         V_usv_y(1 , (ii-1)*s_no+2) = -1;  V_usv_y(1 , (ii-2)*s_no+2) = (1-(1/Vmax));
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+3) = (V_usv_x * x5) ;
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+4) = (V_usv_y * x5) ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+3) = (V_usv_x )' ;          % gradiant
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+4) = (V_usv_y )' ;
         V_usv_x = zeros(1,s_no_mpc); V_usv_y = zeros(1,s_no_mpc); 
      end 
  end
end
  %% Constraint << C1_12 >> and its grad : inequality constraint associated with z should be higer than sea floor (uobstacle-free) 
 % number of inequality constraints for Num_C1_12: auv_n * NumSlot
if C1_12 == 1
          y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + (1 : auv_n * NumSlot)) = (- (Z_AUV_MPC_all * ( A_P0 + B_delt*x5 )) - Z_floor')';
      grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + (1 : auv_n * NumSlot)) = - (Z_AUV_MPC_all * B_delt )';      % gradiant 
end 
  %% Constraint << C5_2 >> and its grad : inequality constraint associated A2U graph 
 % number of inequality constraints for Num_C5_2: 2*NumSlot
 if C5_2 == 1
    sig_n_obj = zeros(3,s_no_mpc); sig_n_obj(1:2,1:2) = -eye(2);  %A2U graph
    for ii = 1 : NumSlot
        for jj = 1: auv_n
            sig_n_obj(1:3,(ii-1)*s_no+((jj-1)*3+(4:6))) = sig_mat(ii,jj)*eye(3);
        end 
        y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + Num_C1_12 + (ii-1)*2+1 ) =  (sig_n_obj * (A_P0 + B_delt*x5) )' * (sig_n_obj * (A_P0 + B_delt*x5) ) - (ds^2 + A2V_tol^2) ;
        y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + Num_C1_12 + (ii-1)*2+2 ) = -(sig_n_obj * (A_P0 + B_delt*x5) )' * (sig_n_obj * (A_P0 + B_delt*x5) ) + (ds^2 - A2V_tol^2) ;
    grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + Num_C1_12 + (ii-1)*2+1 ) =  ((sig_n_obj*B_delt)' * ( (sig_n_obj*A_P0) + (sig_n_obj*B_delt*x5) ) ); 
    grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + Num_C1_12 + (ii-1)*2+2 ) = -((sig_n_obj*B_delt)' * ( (sig_n_obj*A_P0) + (sig_n_obj*B_delt*x5) ) );                          
        sig_n_obj = zeros(3,s_no_mpc); sig_n_obj(1:2,1:2) = -eye(2);
    end  
 end

   %% Constraint << C5_3 >> and its grad : inequality constraint associated A2A graph 
 % number of inequality constraints for Num_C5_3: 2*NumSlot*(auv_n-1)
 if C5_3 == 1
    sig_ij_obj = zeros(3,s_no_mpc);  %A2A graph
    for ii = 1 : NumSlot
        for jj = 1 : auv_n-1
            ij_1 = sig_ij( (jj-1)*2+1 , ii);
            ij_2 = sig_ij( (jj-1)*2+2 , ii);
            sig_ij_obj(1:3, (ii-1)*s_no + (ij_1-1)*3+(4:6)  ) =  eye(3);
            sig_ij_obj(1:3, (ii-1)*s_no + (ij_2-1)*3+(4:6)  ) = -eye(3);
        y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + Num_C1_12 + Num_C5_2 + (ii-1)*2*(auv_n-1)+(jj-1)*2+1 ) =  (sig_ij_obj * (A_P0 + B_delt*x5) )' * (sig_ij_obj * (A_P0 + B_delt*x5) ) - (ds^2 + A2V_tol^2) ;
        y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + Num_C1_12 + Num_C5_2 + (ii-1)*2*(auv_n-1)+(jj-1)*2+2 ) = -(sig_ij_obj * (A_P0 + B_delt*x5) )' * (sig_ij_obj * (A_P0 + B_delt*x5) ) + (ds^2 - A2V_tol^2) ;
    grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + Num_C1_12 + Num_C5_2 + (ii-1)*2*(auv_n-1)+(jj-1)*2+1 ) =  ((sig_ij_obj*B_delt)' * ( (sig_ij_obj*A_P0) + (sig_ij_obj*B_delt*x5) ) ); 
    grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + Num_C1_12 + Num_C5_2 + (ii-1)*2*(auv_n-1)+(jj-1)*2+2 ) = -((sig_ij_obj*B_delt)' * ( (sig_ij_obj*A_P0) + (sig_ij_obj*B_delt*x5) ) );                          
        sig_ij_obj = zeros(3,s_no_mpc);
        end
    end  
 end

  %% Equality Constraint 
yeq = [];
gradyeq = [];
