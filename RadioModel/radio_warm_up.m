function lattice = radio_warm_up(lattice,Kernel,warm_epoch,disp_sign,Beta)
    edge_size = length(lattice);
    Potential = conv2(lattice,Kernel,'same');
    for warm_index = 1:1:warm_epoch
        if mod(warm_index,warm_epoch/100) == 0 && disp_sign == 1
            fprintf("Warm Ratio = %f\n",warm_index/warm_epoch);
        end
        x_change = randi([1,edge_size]);
        y_change = randi([1,edge_size]);
        energy_diff = 2*Potential(x_change,y_change)*(-2*lattice(x_change,y_change));%factor 2 means the energy it give to others
        p_change = exp(-Beta*energy_diff);
        if rand < p_change
            Potential = Potential - 2*lattice(x_change,y_change)* ...
            Kernel(edge_size+2-x_change:2*edge_size+1-x_change,edge_size+2-y_change:2*edge_size+1-y_change);
            lattice(x_change,y_change) = -lattice(x_change,y_change);
        end
    end
end

