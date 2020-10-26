% The Queen in title means she can get many friends even in a distant.

%% Set up parametres
edge_size = 50;
J = 1;
K = 10;
Beta = 0.02;
warm_epoch = 10^2*edge_size^2;
warm_epoch_inside = 1;
mc_epoch = 1e4;
mc_result = zeros(mc_epoch,1);
accept_count = 0;

%% Initialization
lattice = 2*double(rand(edge_size)<0.5)-1;

%% Warming
lattice = warm_up(lattice,J,K,Beta,edge_size,warm_epoch,1);

%% Measuring
for mc_index = 1:1:mc_epoch
    if mod(mc_index,mc_epoch/100)==0
        fprintf("MC_ratio=%f\n",mc_index/mc_epoch);
    end
    x_try = randi([1,edge_size]);
    y_try = randi([1,edge_size]);
    lattice_old = lattice;
    lattice_new = lattice;
    lattice_new(x_try,y_try) = - lattice_new(x_try,y_try);
    [J_factor_new,C_factor_new,energy_new] = measure_C_factor_and_energy(lattice_new,connection,J,K);
    [J_factor_old,C_factor_old,energy_old] = measure_C_factor_and_energy(lattice_old,connection,J,K);
    accept_ratio = exp(-Beta*((energy_new-energy_old)));
    if rand<accept_ratio
        lattice = lattice_new;
        accept_count = accept_count + 1;
    end
    mc_result(mc_index) = sum(lattice,'all');
end
    
%plot(mc_result);
fprintf("Accept_Ratio=%f\n",accept_count/mc_epoch);
sample = mc_result - mean(mc_result);
%% Draw the correlatio curve
sample_len = length(sample);
corr_point = 100;
step = mc_epoch/(10*corr_point);
c02 = mean(sample.^2);
corr_data = zeros(corr_point,1);
for corr_index = 1:1:corr_point
    corr_len = corr_index * step;
    corr_temp = (sample(1:sample_len-corr_len)-mean(sample(1:sample_len-corr_len))).*(sample(1+corr_len:sample_len)-mean(sample(1+corr_len:sample_len)));
    corr_data(corr_index) = mean(corr_temp)/c02;
end
plot(step:step:corr_point*step,corr_data,'blue')
hold on
