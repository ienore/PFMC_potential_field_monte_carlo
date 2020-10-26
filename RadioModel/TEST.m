%% Set up Parameters
edge_size = 50;
A = 0.01;
alpha = 1;
A_c = 0.1;
alpha_c = 2;
Beta = 0.02;
warm_epoch = 100*edge_size.^2;
max_friend = 10;
mc_epoch = 1000;
    mc_energy = zeros([mc_epoch,1]);
    mc_J_factor = zeros([mc_epoch,1]);
    mc_C_factor = zeros([mc_epoch,1]);
    mc_diff = zeros([mc_epoch,1]);
%% Initialize Working Matrixs
lattice = 2*double(rand(edge_size)<0.5)-1;
accept_count = 0;
x = (1:1:(2*edge_size+1))-(edge_size+1);
y = (1:1:(2*edge_size+1))-(edge_size+1);
[X,Y] = meshgrid(x,y);
R = sqrt((X.^2 + Y.^2));
R((edge_size+1),(edge_size+1))=1;
Kernel = -A.* R.^(-alpha);
Kernel((edge_size+1),(edge_size+1))=0;
%% Warming
lattice = radio_warm_up(lattice,Kernel,warm_epoch,1,Beta);
%% Set up Connections
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
        pos_mat = rand(edge_size)< (A_c./sqrt((dis_mat_x - x_index).^2 + (dis_mat_y - y_index).^2).^(alpha_c));
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
%% MC mearsuring
Potential = conv2(lattice,Kernel,'same');
energy = sum(Potential.*lattice,'all');
for mc_index = 1:1:mc_epoch
    if mod(mc_index,mc_epoch/100) == 0 
        fprintf("MC Ratio = %f\n",mc_index/mc_epoch);
    end
    x_change = randi([1,edge_size]);
    y_change = randi([1,edge_size]);
    lattice_new = lattice;
    lattice_new(x_change,y_change) = -lattice_new(x_change,y_change);
    energy_old = energy;
    energy_new = energy_radio_model(lattice_new,Kernel);
    energy_diff = energy_new - energy_old;
    p_change = exp(-Beta*energy_diff);
    if rand < p_change
        lattice(x_change,y_change) = -lattice(x_change,y_change);
        energy = energy_new;
    else
        energy = energy_old;
    end
    mc_energy(mc_index) = energy;
    [mc_C_factor(mc_index),mc_J_factor(mc_index)] = measure_C_factor_and_J_factor(lattice,connection);
    mc_diff(mc_index) = energy - energy_radio_model(lattice,Kernel);
end
%% Plot the result
mc_energy_norm = (mc_energy - mean(mc_energy))/sqrt(var(mc_energy));
mc_C_factor_norm = (mc_C_factor - mean(mc_C_factor))/sqrt(var(mc_C_factor));
mc_J_factor_norm = (mc_J_factor - mean(mc_J_factor))/sqrt(var(mc_J_factor));
scatter(mc_J_factor_norm,mc_energy_norm,10,'filled','blue')
hold on
scatter(mc_C_factor_norm,mc_energy_norm,10,'filled','red')
