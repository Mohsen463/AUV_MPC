%% P4: A2A Adjacency matrix graph optimization    
ii = auv_n; a2a_num = 0;
while ii > 1
    ii  = ii-1;
     a2a_num = a2a_num + ii;   % the a2a_num is a 1_by_(the number of possible A2A edges without repeatation) matrix and edges are assigned to entries of this matrix.
end
A2A_d = zeros(3,a2a_num,NumSlot);
for ii = 1 : NumSlot 
    zz = 1;
    for jj = 1 : auv_n - 1
        for kk = jj + 1 : auv_n
            A2A_d (2,zz,ii) = jj;
            A2A_d (3,zz,ii) = kk;
            A2A_d (1,zz,ii) = norm(Trj_MPC(:,ii+1,jj+1)'-Trj_MPC(:,ii+1,kk+1)');
            zz = zz+1;
        end
    end
end
tic
sig_ij_mat = zeros(NumSlot,a2a_num);
for ii = 1 : NumSlot
    fun4 = @(x4)quadobj4(x4,ii,A2A_d);
    nonlconstr4 = @(x4)quadconstr4(x4,auv_n,a2a_num);
    x0_4 = zeros(a2a_num,1); % Column vector
    [x4,fval,eflag,output,lambda] = fmincon(fun4,x0_4,[],[],[],[],[],[],nonlconstr4);
    sig_ij_mat(ii,:) = x4';
end
toc
sig_ij = zeros(2*(auv_n-1),NumSlot);
    for ii = 1 : NumSlot
        zzz = 0;
        for jj = 1 : a2a_num  
            if abs( sig_ij_mat(ii,jj) - 1 ) < 0.05
               sig_ij(1+zzz,ii) = A2A_d (2,jj,ii);
               sig_ij(2+zzz,ii) = A2A_d (3,jj,ii);
               zzz = zzz+2 ;
            end  
        end  
    end 
    
figure(4),hold on 
for ii = 1 : NumSlot
    for kk = 1 : auv_n-1
        jj_1 = sig_ij((kk-1)*2+1,ii) + 1;
        jj_2 = sig_ij((kk-1)*2+2,ii) + 1;
    % scatter3(Trj_MPC(1,ii+1,jj),Trj_MPC(2,ii+1,jj),Trj_MPC(3,ii+1,jj),15,"m","filled"), hold on   
        plot3([Trj_MPC(1,ii+1,jj_1), Trj_MPC(1,ii+1,jj_2)],[Trj_MPC(2,ii+1,jj_1),Trj_MPC(2,ii+1,jj_2)], [Trj_MPC(3,ii+1,jj_1),Trj_MPC(3,ii+1,jj_2)],"LineWidth",0.5,"Color","y"), hold on
    end 
end
