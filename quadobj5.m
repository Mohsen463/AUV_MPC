function [y,grady] = quadobj5(x5,delta,P_0,W,usv_final,A_g_MPC_all,B_g_MPC_all,Z_AUV_MPC,P_usv_MPC)
         y = W(1,1) * (x5'*x5) ...  %  input energy effort
           + W(1,2) * Z_AUV_MPC*( A_g_MPC_all*P_0 + delta*B_g_MPC_all*x5 ) ...  % deep dive - minimizing distance to sea floor
           + W(1,3) * ( (P_usv_MPC*( A_g_MPC_all*P_0 + delta*B_g_MPC_all*x5)) - usv_final )' * ( ( P_usv_MPC*( A_g_MPC_all*P_0 + delta*B_g_MPC_all*x5) - usv_final ) ) ;  %  minimizing final USV path to a given point
         % - W(1,3) * ( P_usv_MPC*( A_g_MPC_all*P_0 + delta*B_g_MPC_all*x))' * ( P_usv_MPC*( A_g_MPC_all*P_0 + delta*B_g_MPC_all*x) ) ...  %  maximizing the path for more exploration in the given time

  %% gradiant of the objective function (just take derivitave with respect to x)
if nargout > 1
     grady = 2*W(1,1)*x5 ...
           + W(1,2) * (Z_AUV_MPC * (delta*B_g_MPC_all))' ...
           + W(1,3) * 2 * (  ((P_usv_MPC)*(delta*B_g_MPC_all))'  *  (  ((P_usv_MPC)*(A_g_MPC_all*P_0)) + ((P_usv_MPC)*(delta*B_g_MPC_all*x5)) - usv_final )  ) ; %  maximizing the path for more exploration in the given time
         % - W(1,3) * 2 * (  ((P_usv_MPC)*(delta*B_g_MPC_all))'  *  (  ((P_usv_MPC)*(A_g_MPC_all*P_0)) + ((P_usv_MPC)*(delta*B_g_MPC_all*x))  )  )...  %  maximizing the path for more exploration in the given time

end