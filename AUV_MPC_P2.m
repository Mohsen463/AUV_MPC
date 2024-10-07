%%%%%%%%%%%            In the name of GoD    % %%%%%%%%%
 clc, clear,  close all
% Creat Map
depth = 400;  OceanEnvironment;
[X_max, Y_max] = size(Z);
%% Inputs: Path length, sample time, parameters and hyper_parameters
auv_n = 5;  NumSlot = 20; % UAV numbers and number of slots 
  Vmax = 1.5;  delta = 100 ;  ds = 150;  % speed limit, sample time, and A2A distance
  usv_head = pi/4; usv_final = [Y_max-50;X_max-50;0]; W = [0.1 250 2 ] ;  A2V_tol = ds/5 ; 
%% states
s_no = (1+auv_n)*3; s_no_mpc = s_no * NumSlot;
Z_floor = depth * ones (1,auv_n*NumSlot);
%% state space matrixes for MPC modeling 
 a_auv = eye(3);     a_usv = eye(3); a_usv(3,3)=0;    A_auv = eye(3*auv_n );    A_g = blkdiag(a_usv,A_auv);  % matrix A of the pack
 A_g_MPC = A_g^(NumSlot) ;   % matrix A of the pack in all steps gives the last state
 A_g_MPC_all = zeros(s_no_mpc,s_no);    for ii = 1:NumSlot, A_g_MPC_all( (ii-1)*s_no + 1 : (ii-1)*s_no + s_no , 1 : s_no) = A_g^(ii); end    % matrix A of the pack in all steps gives all states

 b_auv = eye(3);     b_usv = eye(3); b_usv(3,3)=0;    B_auv = eye(3*auv_n );    B_g = blkdiag(b_usv,B_auv);  % matrix B of the pack
 B_g_MPC = zeros(s_no_mpc,s_no);    for ii = 1:NumSlot, B_g_MPC( (ii-1)*s_no + 1 : (ii-1)*s_no + s_no , 1 : s_no ) = A_g^(NumSlot-ii) * B_g; end    % matrix B of the pack in all steps gives the last state
 B_g_MPC_all = zeros(s_no_mpc,s_no_mpc);  for ii = 1:NumSlot 
                                              for jj = 1:ii
                                                  B_g_MPC_all( (ii-1)*s_no + 1 : (ii-1)*s_no + s_no , (jj-1)*s_no + 1 : (jj-1)*s_no + s_no ) = A_g^(ii-jj) * B_g; 
                                              end 
                                          end     % matrix B of the pack in all steps gives all states 

 z_auv = [0 0 1];    Z_auv = zeros(1,s_no);   for ii = 1:auv_n, Z_auv(1,3+(ii-1)*3+3) = 1; end   % extract z's of the pack, e.g., for Usv+2 auv : [ 0 0, 0 0 1, 0 0 1]_1by9 
 Z_AUV_MPC = zeros(1,s_no_mpc);  for ii = 1:NumSlot, Z_AUV_MPC( 1 , (ii-1)*s_no + 1 : (ii-1)*s_no + s_no) = Z_auv; end % extract z's of the pack in all steps, e.g., for Usv+2 auv and 3 slots: [ 0 0, 0 0 1, 0 0 1, 0 0, 0 0 1, 0 0 1, 0 0, 0 0 1, 0 0 1 ]_1 by (3*9)
 Z_AUV_MPC_all = zeros(auv_n*NumSlot,s_no_mpc);  for ii = 1:NumSlot, for jj = 1:auv_n, Z_AUV_MPC_all( (ii-1)*auv_n + jj , (ii-1)*s_no + 3 + (jj-1)*3+3) = 1; end, end % extract z's of the pack in all steps, e.g., for Usv+2 auv and 3 slots: [ 0 0, 0 0 1, 0 0 0, 0 0, 0 0 0, 0 0 0, 0 0, 0 0 0, 0 0 0;
        % see -->                                                                                                                                                                                                                                       % 0 0, 0 0 0, 0 0 1, 0 0, 0 0 0, 0 0 0, 0 0, 0 0 0, 0 0 0... ]_6 by (3*9)

 p_usv = zeros(3,s_no); p_usv(1,1) = 1; p_usv(2,2) = 1;  % extract P's of the SUV in the pack
 P_usv_MPC = zeros(3,s_no_mpc);       P_usv_MPC( 1:3 , (NumSlot-1)*s_no + 1 : (NumSlot-1)*s_no + s_no) = p_usv;      % extract P's of the SUV in the last step
 P_usv_MPC_all = zeros(3,s_no_mpc);   for ii = 1:NumSlot, P_usv_MPC_all( 1:3 , (ii-1)*s_no + 1 : (ii-1)*s_no + s_no) = p_usv; end   % extract P's of the SUV in all steps

 p_auv = zeros(3*auv_n,s_no); for ii = 1:auv_n, p_auv((ii-1)*3+1:(ii-1)*3+3,3+(ii-1)*3+1:(ii-1)*3+6)=eye(3); end  % extract P's of the AUVs in the pack
 P_auv_MPC = zeros(3*auv_n,s_no_mpc);  P_auv_MPC( 1:3*auv_n , (NumSlot-1)*s_no + 1 : (NumSlot-1)*s_no + s_no) = p_auv;      % extract P's of the AUVs in the last step
 P_auv_MPC_all = zeros(3*auv_n*NumSlot,s_no_mpc);   for ii = 1:NumSlot, P_auv_MPC_all( (ii-1)*(3*auv_n)+1:(ii-1)*(3*auv_n)+3*auv_n  , (ii-1)*s_no + 1 : (ii-1)*s_no + s_no) = p_auv; end   % extract P's of the AUVs in all steps
 
 P_pack_MPC = zeros(3*(auv_n+1),s_no_mpc);   P_pack_MPC (1:3,1:s_no_mpc) = P_usv_MPC;   P_pack_MPC (4:3*(auv_n+1),1:s_no_mpc) = P_auv_MPC; 
 %% Start positions
 P_0 = zeros(s_no,1) ; for ii = 1 : 3*auv_n, P_0(ii+3) = (-1)^ii * rand(1)*ds; end 
                        for ii = 1 : auv_n, if P_0(3+(ii-1)*3+3) >= 0, P_0(3+(ii-1)*3+3) = -P_0(3+(ii-1)*3+3); end, end
%% MPC Optimization with relaxed binary variabls: P2
tic 
DONE = 1; % DONE = 0 for updating Z_floor, DONE = 1 for skipping Z_floor
while DONE < 2
% Linear or Quadratic Objective with Quadratic Constraints
fun = @(x)quadobj(x,P_0,W,usv_final,A_g_MPC_all,B_g_MPC_all,Z_AUV_MPC,P_usv_MPC,delta);
nonlconstr = @(x)quadconstr(x,ds,Vmax,P_0,Z_floor,usv_head,auv_n,s_no,NumSlot,s_no_mpc,A_g_MPC_all,B_g_MPC_all,Z_AUV_MPC_all,P_usv_MPC,delta);
x0 = zeros(s_no_mpc,1); % Column vector
[x,fval,eflag,output,lambda] = fmincon(fun,x0,[],[],[],[],[],[],nonlconstr);
x1=x; 
  for ii = 1 : s_no_mpc, if x(ii,1) > Vmax, x(ii,1) = Vmax; end, end
toc
%% Trajectory based on determied input jerks

P_traj = zeros(3,s_no_mpc);
Trj_MPC = zeros(3,NumSlot+1,auv_n+1);
for jj = 1 : auv_n, Trj_MPC(:,1,jj+1) = P_0(3+(jj-1)*3+1:3+(jj-1)*3+3);  end
for ii = 1 : NumSlot
    for jj = 1 : auv_n + 1   % the jj_th vehicle (USV or AUVs)
        P_traj(1:3 , (ii-1)*s_no+(jj-1)*(3)+1 : (ii-1)*s_no+(jj-1)*3+3) = eye(3);
        Trj_MPC(1:3,ii+1,jj) = P_traj * ( A_g_MPC_all*P_0 + delta*B_g_MPC_all*x) ;
        P_traj = zeros(3,s_no_mpc);
    end  
end
%% Plotting
if DONE == 0
   figure(3), mesh(Z), zlim([-depth 0]), hold on 
  for ii = 1 : NumSlot + 1
     for jj = 1 : auv_n + 1
        figure(3)
        if ii == 1
           if jj == 1
              scatter3(Trj_MPC(1,ii,jj),Trj_MPC(2,ii,jj),Trj_MPC(3,ii,jj),20,"b","filled"), hold on
           else  
              scatter3(Trj_MPC(1,ii,jj),Trj_MPC(2,ii,jj),Trj_MPC(3,ii,jj),15,"r","filled"), hold on
           end  
        else 
           % CoLoR(uu,:) =  rand(1,3);
           if jj == 1
              scatter3(Trj_MPC(1,ii,jj),Trj_MPC(2,ii,jj),Trj_MPC(3,ii,jj),7.5,"b","filled"), hold on   
              plot3(Trj_MPC(1,ii-1:ii,jj),Trj_MPC(2,ii-1:ii,jj),Trj_MPC(3,ii-1:ii,jj),"LineWidth",1.2,"Color","c"), hold on
           else  
              scatter3(Trj_MPC(1,ii,jj),Trj_MPC(2,ii,jj),Trj_MPC(3,ii,jj),5,"r","filled"), hold on   
              plot3(Trj_MPC(1,ii-1:ii,jj),Trj_MPC(2,ii-1:ii,jj),Trj_MPC(3,ii-1:ii,jj),"LineWidth",1,"Color",[.8 .8 .8]), hold on
           end     
        end     
     end    
  end   

  Trj_MPC_1 = Trj_MPC;
  for jj = 2 : auv_n + 1
     for ii = 2 : NumSlot + 1 
        if floor(Trj_MPC(1,ii,jj)) > 0 && floor(Trj_MPC(2,ii,jj)) > 0 && floor(Trj_MPC(1,ii,jj)) < Y_max && floor(Trj_MPC(2,ii,jj)) < X_max
           if Trj_MPC(3,ii,jj) < Z ( floor(Trj_MPC(2,ii,jj)) , floor(Trj_MPC(1,ii,jj)) )
              Z_floor(1,(ii-2)*auv_n+jj-1) = abs (Z ( floor(Trj_MPC(2,ii,jj)) , floor(Trj_MPC(1,ii,jj)) ) + 1);
           end 
        end 
     end   
  end  

end 
DONE = DONE+1
end


Trj_MPC_2 = Trj_MPC;
for jj = 2 : auv_n + 1
    for ii = 2 : NumSlot + 1 
        if floor(Trj_MPC(1,ii,jj)) > 0  &&   floor(Trj_MPC(2,ii,jj)) > 0   &&   floor(Trj_MPC(1,ii,jj)) < Y_max   &&   floor(Trj_MPC(2,ii,jj)) < X_max
           if Trj_MPC(3,ii,jj) < Z ( floor(Trj_MPC(2,ii,jj)) , floor(Trj_MPC(1,ii,jj)) )
              Trj_MPC(3,ii,jj) = Z ( floor(Trj_MPC(2,ii,jj)) , floor(Trj_MPC(1,ii,jj)) ) + 17;
           end 
        end 
    end  
end 

figure(4), mesh(Z), zlim([-depth 0]), hold on 
for ii = 1 : NumSlot + 1
    for jj = 1 : auv_n + 1
        figure(4)
        if ii == 1
           if jj == 1
              scatter3(Trj_MPC(1,ii,jj),Trj_MPC(2,ii,jj),Trj_MPC(3,ii,jj),20,"b","filled"), hold on
           else  
              scatter3(Trj_MPC(1,ii,jj),Trj_MPC(2,ii,jj),Trj_MPC(3,ii,jj),15,"r","filled"), hold on
           end  
        else 
           % CoLoR(uu,:) =  rand(1,3);
           if jj == 1
              scatter3(Trj_MPC(1,ii,jj),Trj_MPC(2,ii,jj),Trj_MPC(3,ii,jj),10,"b","filled"), hold on   
              plot3(Trj_MPC(1,ii-1:ii,jj),Trj_MPC(2,ii-1:ii,jj),Trj_MPC(3,ii-1:ii,jj),"LineWidth",1.5,"Color","c"), hold on
           else  
              scatter3(Trj_MPC(1,ii,jj),Trj_MPC(2,ii,jj),Trj_MPC(3,ii,jj),7.5,"r","filled"), hold on   
              plot3(Trj_MPC(1,ii-1:ii,jj),Trj_MPC(2,ii-1:ii,jj),Trj_MPC(3,ii-1:ii,jj),"LineWidth",1,"Color",[.8 .8 .8]), hold on
           end     
        end     
    end  
end  

%% A2U Adjacency matrix graph optimization
P3
%% A2A Adjacency matrix graph optimization
P4
%% MPC optimization with binary variables known from P3 and P4
P5








