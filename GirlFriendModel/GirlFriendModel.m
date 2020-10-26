% H = - J \sum S_i S_j - K \sum S_i S_j S_k S_l
%Periodic Boundary Condition
%% Set up parametres
edge_size = 50;
J = 1;
K = 0.2;
Beta = 1;
warm_epoch = 1e2 * edge_size^2;
warm_epoch_inside = 10;
mc_epoch = 1 * edge_size^2;
connection_ratio = 1;
%% Initialization
lattice = 2*double(rand(edge_size)<0.5)-1;
energy_old = measure_total_energy(lattice,J,K);
mc_result_energy = zeros([1,mc_epoch]);
mc_result_J_factor = zeros([1,mc_epoch]);
mc_result_C_factor = zeros([1,mc_epoch]);

%% Initialize the connection
for x_index = 1:1:edge_size
    for y_index = 1:1:edge_size
        connection_index = (x_index - 1)*edge_size + y_index;
        x = x_index;
        y = y_index;
        x0 = mod(x-2,edge_size)+1;
        x1 = mod(x,edge_size)+1;
        y0 = mod(y-2,edge_size)+1;
        y1 = mod(y,edge_size)+1;
        choice_list = zeros([4,2]);
        choice_list(1,:) = [x0,y0];
        choice_list(2,:) = [x1,y0];
        choice_list(3,:) = [x0,y1];
        choice_list(4,:) = [x1,y1];
        if  rand < connection_ratio
            choice = choice_list(randi([1,4]),:);
        else
            choice = [0,0];
        end
        connection(connection_index,:) = choice;
    end
end
%% Warming
lattice = warm_up(lattice,J,K,Beta,edge_size,warm_epoch,1);
pcolor(lattice);

%% Mearsuring
for mc_index = 1:1:mc_epoch
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

scatter(mc_result_C_factor + mc_result_J_factor,mc_result_energy)
p = polyfit(mc_result_J_factor,mc_result_energy,1);
corr_mat = corrcoef(mc_result_C_factor+mc_result_J_factor,mc_result_energy);
corr_mat_no_C = corrcoef(mc_result_J_factor,mc_result_energy);
fprintf("\n Slope = %f\n Corrcoef = %f\n WithoutC = %f\n",p(1),corr_mat(1,2),corr_mat_no_C(1,2))