function lattice = warm_up(lattice,J,K,Beta,edge_size,warm_epoch,disp_sign)
    for warm_index = 1:1:warm_epoch
        if mod(warm_index,warm_epoch/100)==0 && disp_sign == 1
            fprintf("Warm Ratio = %f\n",warm_index/warm_epoch);
        end
        x_to_change = randi([1,edge_size]);
        y_to_change = randi([1,edge_size]);
        x = x_to_change;
        y = y_to_change;
        x0 = mod(x-2,edge_size)+1;
        x1 = mod(x,edge_size)+1;
        y0 = mod(y-2,edge_size)+1;
        y1 = mod(y,edge_size)+1;
        energy_j = -J*lattice(x,y)*(lattice(x0,y)+lattice(x1,y)+lattice(x,y0)+lattice(x,y1));
        energy_k = -K*lattice(x,y)*(lattice(x0,y0)+lattice(x0,y1)+lattice(x1,y1)+lattice(x1,y0));
        energy_step =  energy_j + energy_k;
        energy_diff = 2*energy_step;
        if rand < exp(Beta*(energy_diff))
            lattice(x,y) = -lattice(x,y);
        end
    end
end

