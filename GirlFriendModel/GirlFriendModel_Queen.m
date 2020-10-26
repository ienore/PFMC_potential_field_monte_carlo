% The Queen in title means she can get many friends even in a distant.

%% Set up parametres
edge_size = 50;
J = 1;
K = 10;
Beta = 0.05;
warm_epoch = 10^2 * edge_size^2;
warm_epoch_inside = 1;
mc_epoch = 1e4;
mc_result = zeros(mc_epoch,1);

A = 1600;
alpha = 20;
max_friend = 10;
%% Initialization
lattice = 2*double(rand(edge_size)<0.5)-1;
energy_old = measure_total_energy(lattice,J,K);
mc_result_energy = zeros([1,mc_epoch]);
mc_result_J_factor = zeros([1,mc_epoch]);
mc_result_C_factor = zeros([1,mc_epoch]);
connection = zeros([edge_size^2,max_friend*2]);
friend_count = zeros([edge_size^2,1])+1;

%% Initialize the connection
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
    disp(x_index/edge_size);
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
%% Warming
lattice = warm_up(lattice,J,K,Beta,edge_size,warm_epoch,1);
pcolor(lattice);

%% Mearsuring
for mc_index = 1:1:mc_epoch
    mc_result(mc_index) = sum(lattice,'all');
    if mod(mc_index,mc_epoch/100)==0
        fprintf("MC Ratio = %f\n",mc_index/mc_epoch);
        %disp(norm(lattice-1));
    end
    lattice = warm_up(lattice,J,K,Beta,edge_size,warm_epoch_inside,0);
    [mc_result_J_factor(mc_index),mc_result_C_factor(mc_index),mc_result_energy(mc_index)] =  measure_C_factor_and_energy(lattice,connection,J,K);
end

%% Plot the final result
pcolor(lattice)
axis equal

scatter(( mc_result_J_factor)/mean( mc_result_J_factor),mc_result_energy)
hold on
scatter((mc_result_C_factor )/mean(mc_result_C_factor ),mc_result_energy,'r')

p = polyfit(mc_result_J_factor,mc_result_energy,1);
p_c = polyfit(mc_result_C_factor,mc_result_energy,1);
corr_mat = corrcoef(mc_result_C_factor,mc_result_energy);
corr_mat_no_C = corrcoef(mc_result_J_factor,mc_result_energy);
fprintf("\n J_old = %f\t J_eff = %f\n Corrcoef = %f\t WithoutC = %f\n",p(1),p_c(1),corr_mat(1,2),corr_mat_no_C(1,2))
fprintf(" EnhancedRatio = %f\n",corr_mat(1,2)/corr_mat_no_C(1,2))