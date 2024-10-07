%% MPC optimization with binary variables known from P3 and P4 
tic
fun5 = @(x5)quadobj5(x5,delta,P_0,W,usv_final,A_g_MPC_all,B_g_MPC_all,Z_AUV_MPC,P_usv_MPC);
nonlconstr5 = @(x5)quadconstr5(x5,delta,ds,A2V_tol,Vmax,P_0,Z_floor,usv_head,auv_n,s_no,NumSlot,s_no_mpc,A_g_MPC_all,B_g_MPC_all,Z_AUV_MPC_all,P_usv_MPC,sig_mat,sig_ij);
x0_5 = zeros(s_no_mpc,1); % Column vector
[x5,fval,eflag,output,lambda] = fmincon(fun5,x0_5,[],[],[],[],[],[],nonlconstr5);
x5_2 = x5; 
  for ii = 1 : s_no_mpc, if x5(ii,1) > Vmax, x5(ii,1) = Vmax; end, end
toc
%% Trajectory based on determied input jerks

Trj_MPC_P5 = zeros(3,NumSlot+1,auv_n+1);
for jj = 1 : auv_n, Trj_MPC_P5(:,1,jj+1) = P_0(3+(jj-1)*3+1:3+(jj-1)*3+3);  end
P_traj = zeros(3,s_no_mpc);
for ii = 1 : NumSlot
    for jj = 1 : auv_n + 1   % the jj_th vehicle (USV or AUVs)
        P_traj(1:3 , (ii-1)*s_no+(jj-1)*(3)+1 : (ii-1)*s_no+(jj-1)*3+3) = eye(3);
        Trj_MPC_P5(1:3,ii+1,jj) = P_traj * ( A_g_MPC_all*P_0 + delta*B_g_MPC_all*x5) ;
        P_traj = zeros(3,s_no_mpc);
    end  
end
%% Plotting
Trj_MPC_P5_2 = Trj_MPC_P5;
for jj = 2 : auv_n + 1
    for ii = 2 : NumSlot + 1 
        if floor(Trj_MPC_P5(1,ii,jj)) > 0  &&   floor(Trj_MPC_P5(2,ii,jj)) > 0   &&   floor(Trj_MPC_P5(1,ii,jj)) < Y_max   &&   floor(Trj_MPC_P5(2,ii,jj)) < X_max
           if Trj_MPC_P5(3,ii,jj) < Z ( floor(Trj_MPC_P5(2,ii,jj)) , floor(Trj_MPC_P5(1,ii,jj)) )
              Trj_MPC_P5(3,ii,jj) = Z ( floor(Trj_MPC_P5(2,ii,jj)) , floor(Trj_MPC_P5(1,ii,jj)) ) + 17;
           end 
        end 
    end  
end 

   figure(5), mesh(Z), zlim([-depth 0]), hold on 
for ii = 1 : NumSlot + 1
    for jj = 1 : auv_n + 1
        figure(5)
        if ii == 1
           if jj == 1
              scatter3(Trj_MPC_P5(1,ii,jj),Trj_MPC_P5(2,ii,jj),Trj_MPC_P5(3,ii,jj),20,"b","filled"), hold on
           else  
              scatter3(Trj_MPC_P5(1,ii,jj),Trj_MPC_P5(2,ii,jj),Trj_MPC_P5(3,ii,jj),15,"r","filled"), hold on
           end  
        else 
           % CoLoR(uu,:) =  rand(1,3);
           if jj == 1
              scatter3(Trj_MPC_P5(1,ii,jj),Trj_MPC_P5(2,ii,jj),Trj_MPC_P5(3,ii,jj),10,"b","filled"), hold on   
              plot3(Trj_MPC_P5(1,ii-1:ii,jj),Trj_MPC_P5(2,ii-1:ii,jj),Trj_MPC_P5(3,ii-1:ii,jj),"LineWidth",1.5,"Color","c"), hold on
           else  
              scatter3(Trj_MPC_P5(1,ii,jj),Trj_MPC_P5(2,ii,jj),Trj_MPC_P5(3,ii,jj),7.5,"r","filled"), hold on   
              plot3(Trj_MPC_P5(1,ii-1:ii,jj),Trj_MPC_P5(2,ii-1:ii,jj),Trj_MPC_P5(3,ii-1:ii,jj),"LineWidth",1,"Color",[.8 .8 .8]), hold on
           end     
        end     
    end  
end  

% figure(5),hold on 
% for ii = 1 : NumSlot
%     jj = sig(1,ii)+1;
%      scatter3(Trj_MPC_P5(1,ii+1,jj),Trj_MPC_P5(2,ii+1,jj),Trj_MPC_P5(3,ii+1,jj),15,"m","filled"), hold on   
%      plot3([Trj_MPC_P5(1,ii+1,1), Trj_MPC_P5(1,ii+1,jj)],[Trj_MPC_P5(2,ii+1,1),Trj_MPC_P5(2,ii+1,jj)], [Trj_MPC_P5(3,ii+1,1),Trj_MPC_P5(3,ii+1,jj)],"LineWidth",0.5,"Color","c"), hold on
% end 
% 
% figure(5),hold on 
% for ii = 1 : NumSlot
%     for kk = 1 : auv_n-1
%         jj_1 = sig_ij((kk-1)*2+1,ii) + 1;
%         jj_2 = sig_ij((kk-1)*2+2,ii) + 1;
%     % scatter3(Trj_MPC(1,ii+1,jj),Trj_MPC(2,ii+1,jj),Trj_MPC(3,ii+1,jj),15,"m","filled"), hold on   
%         plot3([Trj_MPC_P5(1,ii+1,jj_1), Trj_MPC_P5(1,ii+1,jj_2)],[Trj_MPC_P5(2,ii+1,jj_1),Trj_MPC_P5(2,ii+1,jj_2)], [Trj_MPC_P5(3,ii+1,jj_1),Trj_MPC_P5(3,ii+1,jj_2)],"LineWidth",0.5,"Color","y"), hold on
%     end 
% end