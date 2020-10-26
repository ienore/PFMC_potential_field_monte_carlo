% The Queen in title means she can get many friends even in a distant.

%% Set up parametres
edge_size = 100;
J = -1;
K = 1;
Beta = 0.1;
warm_epoch = 10^2*edge_size^2;
warm_epoch_inside = 10;
mc_epoch = 1e5;
mc_potential = zeros(mc_epoch,1);
mc_N1 = zeros(mc_epoch,1);
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
    accept_ratio = get_accept_ratio(lattice,x_try,y_try,J,K,Beta);
    quasi_potential = log(accept_ratio)/lattice(x_try,y_try);
    N1 = get_N1(lattice,x_try,y_try);
    mc_potential(mc_index) = quasi_potential;
    mc_N1(mc_index) = N1;
    lattice = warm_up(lattice,J,K,Beta,edge_size,warm_epoch_inside,0);
end
scatter(mc_N1,mc_potential,5,'filled','r');