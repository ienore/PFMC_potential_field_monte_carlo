% The Queen in title means she can get many friends even in a distant.

%% Set up parametres
edge_size = 50;
J = 1;
K = 10;
Beta = 0.02;
J_effective = 4;
warm_epoch =  10^2*edge_size^2;
mc_epoch = 1e4;
mc_result = zeros(mc_epoch,1);
accept_count = 0;
max_friend = 10;
%% Warming Using Cluster Algorithm
lattice = 2*double(rand(edge_size)<0.5)-1;
lattice = warm_up(lattice,J,K,Beta,edge_size,warm_epoch,1);

%% Measuring With Cluster Algorithm
for mc_index = 1:1:mc_epoch
    if mod(mc_index,mc_epoch/100)==0
        fprintf("MC_ratio=%f\n",mc_index/mc_epoch);
    end
    if mod(mc_index,1000)==0
        fprintf("Reconnection...\n");
       %% Initialize the connection
        connection = zeros([edge_size^2,max_friend*2]);
        friend_count = zeros([edge_size^2,1])+1;
        dis_mat_x = zeros(edge_size);
        dis_mat_y = zeros(edge_size);
        for x_index = 1:1:edge_size
            for y_index = 1:1:edge_size
                dis_mat_x(x_index,y_index) = x_index;
                dis_mat_y(x_index,y_index) = y_index;
            end
        end
        % Set up  the connection for the first time
        for x_index = 1:1:edge_size
            for y_index = 1:1:edge_size
                connection_index = (x_index - 1)*edge_size + y_index;
                pos_mat = rand(edge_size)< (A./sqrt((dis_mat_x - x_index).^2 + (dis_mat_y - y_index).^2).^(alpha));
                for x_try = 1:1:edge_size
                    if friend_count(connection_index) > max_friend
                       break;
                    end
                    for y_try = 1:1:edge_size
                        if friend_count(connection_index) > max_friend
                            break;
                        end
                        if pos_mat(x_try,y_try) == 1 &&( (x_try ~= x_index) || (y_try ~= y_index))
                            % Add x_try,y_try to x_index,y_index's list
                            friend_index = friend_count(connection_index);
                            connection(connection_index,2*(friend_index-1)+1) = x_try;
                            connection(connection_index,2*(friend_index-1)+2) = y_try;
                            friend_count(connection_index) = friend_count(connection_index)+1;

                            connected_index = (x_try - 1)*edge_size + y_try;
                            friend_index_ed = friend_count(connected_index);
                            connection(connected_index,2*(friend_index_ed-1)+1) = x_index;
                            connection(connected_index,2*(friend_index_ed-1)+2) = y_index;
                            friend_count(connected_index)  = friend_count(connected_index)+1;
                        end
                    end
                end
            end
        end
        % Clean the connection, clear the repeated links
        exam_connection = zeros(edge_size);
        for x_index = 1:1:edge_size
            for y_index = 1:1:edge_size
                exam_connection = zeros(edge_size);
                connected_index = (x_index - 1)*edge_size + y_index;
                for friend_index = 1:1:max_friend
                    x_new = connection(connected_index,2*(friend_index-1)+1);
                    y_new = connection(connected_index,2*(friend_index-1)+2);
                    if x_new ~=0 && y_new ~= 0
                        if exam_connection(x_new,y_new) == 1
                            connection(connected_index,2*(friend_index-1)+1) = 0;
                            connection(connected_index,2*(friend_index-1)+2) = 0;
                        else
                            exam_connection(x_new,y_new) = 1;
                        end
                    end
                end
            end
        end
    end
    x_start = randi(edge_size);
    y_start = randi(edge_size);
    cluster_sign = lattice(x_start,y_start);
    new_include_list = zeros([edge_size^2,2]);
    new_include_list_temp = zeros([edge_size^2,2]);
    new_include_list(1,:) = [x_start,y_start];
    new_list_count = 1;
    new_include_square_temp = zeros(edge_size);
    new_list_count_temp = 0;
    p_cluster = 1 - exp(-2*Beta*J_effective);
    change = zeros(edge_size)+1;
    while(new_list_count > 0)%Propose the change
        for new_index = 1:1:new_list_count% Upgrade the new_include
            x_index = new_include_list(new_index,1);
            y_index = new_include_list(new_index,2);
            change(x_index,y_index) = -1;
            connected_index = (x_index - 1)*edge_size + y_index;
            for friend_index = 1:1:max_friend
                x_friend = connection(connected_index,2*(friend_index-1)+1);
                y_friend = connection(connected_index,2*(friend_index-1)+2);
                if x_friend ~= 0 && y_friend ~= 0 && change(x_friend,y_friend)==1 && new_include_square_temp(x_friend,y_friend)==0 %This friend is not empty
                    if lattice(x_friend,y_friend) == cluster_sign && rand < p_cluster %% Add into the cluster
                        new_list_count_temp = new_list_count_temp + 1;
                        new_include_list_temp(new_list_count_temp,:) = [x_friend,y_friend];
                        new_include_square_temp(x_friend,y_friend) = 1;
                    end
                end
            end
        end
        new_include_list = new_include_list_temp;
        new_list_count = new_list_count_temp;
        new_include_list_temp = zeros(edge_size^2,2);
        new_list_count_temp = 0;
    end
    %Decide whether or not to accept the change:
    lattice_old = lattice;
    lattice_new = lattice.*change;
    [J_factor_new,C_factor_new,energy_new] = measure_C_factor_and_energy(lattice_new,connection,J,K);
    [J_factor_old,C_factor_old,energy_old] = measure_C_factor_and_energy(lattice_old,connection,J,K);
    accept_ratio = exp(-Beta*((energy_new-energy_old)-(C_factor_new-C_factor_old)*(-J_effective) ));
    if rand<accept_ratio
        lattice = lattice_new;
        accept_count = accept_count +1;
    end
    mc_result(mc_index) = sum(lattice,'all');
end
%plot(mc_result)
fprintf("AcceptRatio = %f\n",accept_count/mc_epoch);
sample_A = mc_result - mean(mc_result);
%% Draw the correlatio curve
sample_len = length(sample_A);
corr_point = 100;
step = mc_epoch/(10*corr_point);
c02 = mean(sample_A.^2);
corr_data = zeros(corr_point,1);
for corr_index = 1:1:corr_point
    corr_len = corr_index * step;
    corr_temp = (sample_A(1:sample_len-corr_len)-mean(sample_A(1:sample_len-corr_len))).*(sample_A(1+corr_len:sample_len)-mean(sample_A(1+corr_len:sample_len)));
    corr_data(corr_index) = mean(corr_temp)/c02;
end
plot(step:step:corr_point*step,corr_data,'red')
hold on
