%% P4: Longest tree  - LT
 clc, clear, close all
Nodes_num = 5; 
for ii = 1 : Nodes_num
    Nodes(ii,1) = rand(1) * 100;
    Nodes(ii,2) = rand(1) * 5 + 50;
end
  ii = Nodes_num; a2a_num = 0;
while ii > 1
    ii  = ii-1;
     a2a_num = a2a_num + ii;   % the a2a_num is a 1_by_(the number of possible A2A edges without repeatation) matrix and edges are assigned to entries of this matrix.
end
A2A_d = zeros(3,a2a_num);
    zz = 1;
    for jj = 1 : Nodes_num - 1
        for kk = jj + 1 : Nodes_num
            A2A_d (2,zz) = jj;
            A2A_d (3,zz) = kk;
            A2A_d (1,zz) = norm(Nodes(jj,:)-Nodes(kk,:));
            zz = zz+1;
        end
    end

for iii = 1:2
sig_ij_mat = zeros(1,a2a_num);
    fun_LT = @(xLT)quadobj_LT(xLT,A2A_d,iii);
    nonlconstr_LT = @(xLT)quadconstr_LT(xLT,Nodes_num,a2a_num);
    xLT_0 = zeros(a2a_num,1); % Column vector
    [xLT,fval,eflag,output,lambda] = fmincon(fun_LT,xLT_0,[],[],[],[],[],[],nonlconstr_LT);
    sig_ij_mat(1,:) = xLT';
sig_ij = zeros(2*(Nodes_num-1),1);
        zzz = 0;
        for jj = 1 : a2a_num  
            if abs( sig_ij_mat(1,jj) - 1 ) < 0.05
               sig_ij(1+zzz,1) = A2A_d (2,jj);
               sig_ij(2+zzz,1) = A2A_d (3,jj);
               zzz = zzz+2 ;
            end  
        end  
    
figure(4),hold on 
for ii = 1 : Nodes_num
    scatter(Nodes(ii,1),Nodes(ii,2),20,"k","filled"), hold on 
end

    for kk = 1 : Nodes_num-1
        jj_1 = sig_ij((kk-1)*2+1) ;
        jj_2 = sig_ij((kk-1)*2+2) ;
        if iii == 1 % Longest path
           plot([Nodes(jj_1,1), Nodes(jj_2,1)],[Nodes(jj_1,2), Nodes(jj_2,2)],"LineWidth",0.5,"Color","r"), hold on
        else % Shortest path
           plot([Nodes(jj_1,1), Nodes(jj_2,1)],[Nodes(jj_1,2), Nodes(jj_2,2)],"LineWidth",0.5,"Color","b"), hold on
        end
    end 

end 