edge_size = 100;
J = 1;
K = 0.2;
Beta = 10;
mc_epoch = 10 * edge_size^2;

%% Initialization
lattice = 2*double(rand(edge_size)<0.5)-1;
energy_old = measure_total_energy(lattice,J,K);
mc_result_energy = zeros([1,mc_epoch]);
mc_result_J_factor = zeros([1,mc_epoch]);

%% Start MC
for mc_index = 1:1:mc_epoch
    if mod(mc_index,mc_epoch/100)==0
        fprintf("MC Ratio = %f\n",mc_index/mc_epoch);
        %disp(norm(lattice-1));
    end
    lattice = 2*double(rand(edge_size)<0.5)-1;
    [J_factor,energy] = measure_J_factor_and_energy(lattice,J,K);
    mc_result_energy(mc_index) = energy;
    mc_result_J_factor(mc_index) = J_factor;
end

%% Plot the result
scatter(mc_result_J_factor,mc_result_energy,5,'filled','r')
p = polyfit(mc_result_J_factor,mc_result_energy,1);
fprintf("Slope = %f\n",p(1))