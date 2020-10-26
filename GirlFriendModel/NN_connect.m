%% Set up parametres
edge_size = 50;
J = 1;
K = 0;
Beta = 1;
connection_ratio = 0.1;
mc_epoch = 1 * edge_size^2;

%% Initialization
lattice = 2*double(rand(edge_size)<0.5)-1;
connection = zeros([edge_size^2,2]);
energy_old = measure_total_energy(lattice,J,K);
mc_result_energy = zeros([1,mc_epoch]);
mc_result_C_factor = zeros([1,mc_epoch]);
mc_result_J_factor = zeros([1,mc_epoch]);

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
%% Start MC
for mc_index = 1:1:mc_epoch
    if mod(mc_index,mc_epoch/100)==0
        fprintf("MC Ratio = %f\n",mc_index/mc_epoch);
        %disp(norm(lattice-1));
    end
    lattice = 2*double(rand(edge_size)<0.5)-1;
    [J_factor,C_factor,energy] = measure_C_factor_and_energy(lattice,connection,J,K);
    mc_result_energy(mc_index) = energy;
    mc_result_C_factor(mc_index) = C_factor;
    mc_result_J_factor(mc_index) = J_factor;
end

%% Plot the result
scatter(mc_result_C_factor,mc_result_energy,5)
p = polyfit(mc_result_C_factor,mc_result_energy,1);
corr_mat = corrcoef(mc_result_C_factor,mc_result_energy);
fprintf("\n Slope = %f\n Corrcoef = %f\n",p(1),corr_mat(1,2))
scatter(mc_result_J_factor,mc_result_C_factor,5)



