function [y,yeq,grady,gradyeq] = quadconstr(x,ds,Vmax,P_0,Z_floor,usv_head,auv_n,s_no,NumSlot,s_no_mpc,A_g_MPC_all,B_g_MPC_all,Z_AUV_MPC_all,P_usv_MPC,delta)
% modeling the inequality constraints  like y( 1 , 1:Num_ineq_const)
% modeling the grad of inequality constraints  like y( s_no_mpc , 1:Num_ineq_const) and for each column(constraints) just take derivitave of each constraint in y with respect to x
   C1_4 = 1;          C1_6 = 1;            C1_7 = 1;      C1_14 = 1;       C1_12 = 1;  % if constrant X (i.e., C1_X) is taken into acount it is 1, otherwise 0.
% negative z;   dis(auv(K)-Usv(K))<ds;        speed         Heading          floor

Num_C1_4 =  C1_4 * auv_n*NumSlot;       Num_C1_6 = C1_6 * auv_n;             Num_C1_7 = C1_7 * (auv_n+1)*NumSlot;  % number of inequality constraints each C1_X imposes
Num_C1_14 = C1_14 * 4*NumSlot;            Num_C1_12 = C1_12 * auv_n*NumSlot;     
Num_ineq_const = Num_C1_4  +  Num_C1_6  +  Num_C1_7  +  Num_C1_14  +  Num_C1_12 ;   % number of all inequality constraints  

y = zeros(1, Num_ineq_const);
grady = zeros(s_no_mpc,Num_ineq_const);

A_P0 = A_g_MPC_all*P_0;  B_delt = delta*B_g_MPC_all;
%% Constraint << C1_4 >> and its grad : inequality constraint associated with z should be negative (underwater sea level) 
 % number of inequality constraints Num_C1_4: auv_n * NumSlot
 if C1_4 == 1
        y(  1: auv_n * NumSlot) = (Z_AUV_MPC_all * ( A_P0 + B_delt*x ))'; 
    grady(:,1: auv_n * NumSlot) = (Z_AUV_MPC_all * B_delt )';
 end   

 %% Constraint C1_6 and its grad : unequality constraint associated with distances of AUVs from USV at the final stage be lower than a given distance ds based on the sonar modem/sensor rnge
 % number of inequality constraints Num_C1_6 = auv_n 
 if C1_6 == 1
    P_uav_individual = zeros(3,s_no_mpc);   % C1_6
    for jj = 1 : auv_n
        P_uav_individual(1:3 , (NumSlot-1)*s_no+(jj-1)*3+4 : (NumSlot-1)*s_no+(jj-1)*3+6) = eye(3);
            y(  Num_C1_4 + jj) = 1*  ((P_usv_MPC - P_uav_individual) * ( A_P0 + B_delt*x))' * ((P_usv_MPC - P_uav_individual) * ( A_P0 + B_delt*x)) - ds^2 ;
        grady(:,Num_C1_4 + jj) = 2*( ((P_usv_MPC - P_uav_individual)*B_delt)'  *  (  ((P_usv_MPC - P_uav_individual)*A_P0) + ((P_usv_MPC - P_uav_individual)*(B_delt*x))  )  );
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
              y(  Num_C1_4 + Num_C1_6 + (ii-1)*(auv_n + 1)+jj ) =  (V_individual * x)' * (V_individual * x) - Vmax^2 ;        %%% NumSlot must be (auv_n + 1)
          grady(:,Num_C1_4 + Num_C1_6 + (ii-1)*(auv_n + 1)+jj ) =  (V_individual)'     * (V_individual * x) ;  % gradiant     %%% NumSlot must be (auv_n + 1)
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
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+1 ) = (V_usv_x * x)  -  (Vmax*cos(usv_head)) ; % for the first step time specify max of input speed in x & y direction based on the initial heading 
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+2 ) = (V_usv_y * x)  -  (Vmax*sin(usv_head)) ; % for the first step time specify min of input speed in x & y direction be zero (to be positive) to keep in the first quarter of cartesian coordinates 
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+3 ) = -(V_usv_x * x) +  (Vmax*cos(1.5*usv_head)) ; 
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+4 ) = -(V_usv_y * x) +  (Vmax*sin(usv_head/2)) ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+1) =  (V_usv_x )'  ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+2) =  (V_usv_y )'  ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+3) = -(V_usv_x )'  ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+4) = -(V_usv_y )'  ;
         V_usv_x = zeros(1,s_no_mpc); V_usv_y = zeros(1,s_no_mpc); 
      else
         V_usv_x(1 , (ii-1)*s_no+1) = 1;  V_usv_x(1 , (ii-2)*s_no+1) = -(1+(1/Vmax)); % for modeling Vx[k] - Vx[k-1] < Vx[k-1]/Vmax
         V_usv_y(1 , (ii-1)*s_no+2) = 1;  V_usv_y(1 , (ii-2)*s_no+2) = -(1+(1/Vmax));
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+1) =  (V_usv_x * x) ;
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+2) =  (V_usv_y * x) ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+1) =  (V_usv_x )' ;        % gradiant
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+2) =  (V_usv_y )' ;
         V_usv_x = zeros(1,s_no_mpc); V_usv_y = zeros(1,s_no_mpc);

         V_usv_x(1 , (ii-1)*s_no+1) = -1;  V_usv_x(1 , (ii-2)*s_no+1) = (1-(1/Vmax));  % for modeling Vx[k] - Vx[k-1] > -Vx[k-1]/Vmax
         V_usv_y(1 , (ii-1)*s_no+2) = -1;  V_usv_y(1 , (ii-2)*s_no+2) = (1-(1/Vmax));
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+3) = (V_usv_x * x) ;
             y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+4) = (V_usv_y * x) ;
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+3) = (V_usv_x )' ;          % gradiant
         grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + (ii-1)*4+4) = (V_usv_y )' ;
         V_usv_x = zeros(1,s_no_mpc); V_usv_y = zeros(1,s_no_mpc); 
      end 
  end
end
  %% Constraint << C1_12 >> and its grad : inequality constraint associated with z should be higer than sea floor (uobstacle-free) 
 % number of inequality constraints for Num_C1_12: auv_n * NumSlot
if C1_12 == 1
          y(  Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + (1 : auv_n * NumSlot)) = (- (Z_AUV_MPC_all * ( A_P0 + B_delt*x )) - Z_floor')';
      grady(:,Num_C1_4 + Num_C1_6 + Num_C1_7 + Num_C1_14 + (1 : auv_n * NumSlot)) = - (Z_AUV_MPC_all * B_delt )';      % gradiant 
end 
yeq = [];
gradyeq = [];
