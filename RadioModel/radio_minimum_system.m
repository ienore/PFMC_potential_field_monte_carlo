%% Set up Parameters
edge_size = 50;
A = -0.01;
alpha = 0.2;
gamma = 0;
Beta = 0.02;
warm_epoch = 100*edge_size.^2;

%% Initialize Working Matrixs
lattice = 2*double(rand(edge_size)<0.5)-1;
accept_count = 0;
x = (1:1:(2*edge_size+1))-(edge_size+1);
y = (1:1:(2*edge_size+1))-(edge_size+1);
[X,Y] = meshgrid(x,y);
R = sqrt((X.^2 + Y.^2));
R((edge_size+1),(edge_size+1))=1;
Kernel = A.* R.^(-alpha).*exp(-Beta*R);
Kernel((edge_size+1),(edge_size+1))=0;
%% Warming
lattice = radio_warm_up(lattice,Kernel,warm_epoch,1,Beta);


pcolor(lattice);
