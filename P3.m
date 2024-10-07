%% P3: A2U Adjacency matrix graph optimization    
USV_AUV_d = zeros(NumSlot,auv_n);
for ii = 1 : NumSlot 
    for jj = 1 : auv_n 
         USV_AUV_d (ii,jj) = norm(Trj_MPC(:,ii+1,1)'-Trj_MPC(:,ii+1,jj+1)');
    end
end
 tic
for ii = 1 : NumSlot
    fun3 = @(x3)quadobj3(x3,ii,USV_AUV_d);
    nonlconstr3 = @(x3)quadconstr3(x3,auv_n);
    x0_3 = zeros(auv_n,1); % Column vector
    [x3,fval,eflag,output,lambda] = fmincon(fun3,x0_3,[],[],[],[],[],[],nonlconstr3);
    sig_mat(ii,:) = x3';
    for jj = 1 : auv_n
        if abs( sig_mat(ii,jj) - 1 ) < 0.05
            sig(1,ii) = jj;
        end
    end
end
 toc
figure(4),hold on 
for ii = 1 : NumSlot
    jj = sig(1,ii)+1;
     scatter3(Trj_MPC(1,ii+1,jj),Trj_MPC(2,ii+1,jj),Trj_MPC(3,ii+1,jj),15,"m","filled"), hold on   
     plot3([Trj_MPC(1,ii+1,1), Trj_MPC(1,ii+1,jj)],[Trj_MPC(2,ii+1,1),Trj_MPC(2,ii+1,jj)], [Trj_MPC(3,ii+1,1),Trj_MPC(3,ii+1,jj)],"LineWidth",0.5,"Color","c"), hold on
end 
