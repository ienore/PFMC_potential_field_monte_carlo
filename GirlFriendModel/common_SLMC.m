% H = - J \sum S_i S_j - K \sum S_i S_j S_k S_l
%Periodic Boundary Condition
%% Set up parametres
edge_size = 100;
J = 0;
K = 5;
Beta = 1;
warm_epoch = 1e2 * edge_size^2;
warm_epoch_inside = 10;
mc_epoch = 10 * edge_size^2;

%% Initialization
lattice = 2*double(rand(edge_size)<0.5)-1;
energy_old = measure_total_energy(lattice,J,K);
mc_result_energy = zeros([1,mc_epoch]);
mc_result_J_factor = zeros([1,mc_epoch]);
%% Warming
lattice = warm_up(lattice,J,K,Beta,edge_size,warm_epoch,1);


%% Mearsuring
for mc_index = 1:1:mc_epoch
    if mod(mc_index,mc_epoch/100)==0
        fprintf("MC Ratio = %f\n",mc_index/mc_epoch);
        %disp(norm(lattice-1));
    end
    lattice = warm_up(lattice,J,K,Beta,edge_size,warm_epoch_inside,0);
    [mc_result_J_factor(mc_index),mc_result_energy(mc_index)] =  measure_J_factor_and_energy(lattice,J,K);
end

%% Plot the final result
pcolor(lattice)
axis equal

scatter(mc_result_J_factor,mc_result_energy)
p = polyfit(mc_result_J_factor,mc_result_energy,1);
fprintf("Slope = %f\n",p(1))