%% Set up Parameters
zjy_index = 1;
T_hop = 1.0;
NumInEdge = 8;
NumOfVertexs = NumInEdge^2;
K = Get_K(NumInEdge);

Uene = 2;
Miu = Uene/2;
Beta = 2;
D_Tau = 0.2;
TempSlice = Beta/D_Tau;
lambda = 2.0*atanh(sqrt(tanh(D_Tau*Uene/4.0)));
NumOfWarm = 100;
%NumOfWarm = 10;
NumOfEpoch = 100;
Sigma = double(rand([TempSlice,NumOfVertexs])>0.5)*2.0-1.0;%RandomInit
N_wrap = 10;
N_cut = 5.0;
id_mat = eye(NumOfVertexs);
%% Prepare Matrixs
mea_x = 1;
mea_y = 1;
mea_result = zeros([1,TempSlice*NumOfEpoch*2]);
mea_result_auxi = zeros([1,(TempSlice-1)*NumOfEpoch*2]);
%mc_potential = zeros([1,(TempSlice-1)*NumOfEpoch*2*NumOfVertexs]);
%mc_NN = zeros([1,(TempSlice-1)*NumOfEpoch*2*NumOfVertexs]);
MAX_NN  = 4;
mc_potential = zeros([1,NumOfEpoch*NumOfVertexs]);
mc_NN = zeros([NumOfEpoch*NumOfVertexs,TempSlice*MAX_NN]);
count_list = zeros([1,TempSlice])+1;
count = 1.0;
count_new = 1;

%% Warming up
Sigma = WarmUp(zjy_index,N_wrap,Sigma,id_mat,NumInEdge,NumOfWarm,NumOfEpoch,K,TempSlice,NumOfVertexs,Miu,Uene,D_Tau,lambda,T_hop);

%% MC Measuring
for epoch_index = 1:1:NumOfEpoch
    if mod(zjy_index,8) == 1 && mod(epoch_index,NumOfEpoch/1000)==0
        fprintf("D_Tau = %f,MC_Ratio = %f\n",D_Tau,epoch_index/NumOfEpoch);
    end
    for reverse_sign = 0:1:1
    for time_index_mother = 2:1:TempSlice
        if reverse_sign == 0
            time_index = time_index_mother;
            if mod(time_index_mother,N_wrap) == 1 || time_index_mother == 2
                green_L_up = Get_G_L(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                green_L_down = Get_G_L(-1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
            B_trans_up = Get_B_L2(1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_up = Get_B_L2_inv(1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_down = Get_B_L2(-1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_down = Get_B_L2_inv(-1,time_index,time_index-1,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            green_L_up = B_trans_up*green_L_up*B_trans_inv_up;
            green_L_down = B_trans_down*green_L_down*B_trans_inv_down;
            end
        else
            time_index = TempSlice - time_index_mother + 1;
            if mod(time_index_mother,N_wrap) == 1
                green_L_up = Get_G_L(1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
                green_L_down = Get_G_L(-1.0,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            else
            B_trans_up = Get_B_L2(1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_up = Get_B_L2_inv(1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_down = Get_B_L2(-1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            B_trans_inv_down = Get_B_L2_inv(-1,time_index+1,time_index,NumOfVertexs,Sigma,D_Tau,lambda,TempSlice,K,T_hop,Miu,Uene);
            green_L_up = B_trans_inv_up*green_L_up*B_trans_up;
            green_L_down = B_trans_inv_down*green_L_down*B_trans_down;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEASURE (Run :(TempSlice-1)*NumOfEpoch*2) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       	
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MEASURE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for site_index = 1:1:NumOfVertexs
            delta_up = zeros(NumOfVertexs);
            delta_up(site_index,site_index) = exp(-2*lambda*Sigma(time_index,site_index))-1;
            delta_down = zeros(NumOfVertexs);
            delta_down(site_index,site_index) = exp(2*lambda*Sigma(time_index,site_index))-1;
            R_up = 1 + delta_up(site_index,site_index)*(1-green_L_up(site_index,site_index));
            R_down = 1 + delta_down(site_index,site_index)*(1-green_L_down(site_index,site_index));
            PosToChange = R_up*R_down;
            %%% Load Data
            if time_index == 1
                mc_potential(count) = log(PosToChange)/Sigma(time_index,site_index);
                for m_index = 1:1:MAX_NN
                    for time_NN_index = 1:1:TempSlice
                        mc_NN(count,(m_index-1)*TempSlice+time_NN_index) = get_m_th_NN(Sigma,NumInEdge,site_index,time_NN_index,m_index);
                    end
                end
                count = count +1;
            end
            %%% Finished Data
            if rand < PosToChange
                Sigma(time_index,site_index) = -Sigma(time_index,site_index);
                green_L_up = green_L_up - 1.0/R_up * green_L_up * delta_up * (id_mat - green_L_up);
                green_L_down = green_L_down - 1.0/R_down * green_L_down * delta_down * (id_mat - green_L_down);
            end
        end
    end
    end
end



%% Deal with the data
factor_list = rand([TempSlice*MAX_NN,1]);
N_try = 1e5;
step = 0.1;
r2_old = 0;
for try_index = 1:1:N_try
    if mod(try_index,N_try/100)==0
        fprintf("Try_ratio = %f\n",try_index/N_try);
    end
    site_change = randi([1,TempSlice*MAX_NN]);
    factor_list_new = factor_list;
    factor_list_new(site_change) = factor_list_new(site_change) + step * (rand-0.5);
    new_factor = mc_NN * factor_list_new;
    corr_mat = corrcoef(mc_potential,new_factor);
    r2 = abs(corr_mat(1,2)/sqrt(corr_mat(1,1)*corr_mat(2,2)));
    if r2 > r2_old
        factor_list = factor_list_new;
       	disp(r2);
        r2_old = r2;
    end
end
scatter(new_factor,mc_potential)
%% Plot the final result





