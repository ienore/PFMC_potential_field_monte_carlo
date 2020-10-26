function ratio = get_accept_ratio(lattice,x,y,J,K,Beta)
    edge_size = length(lattice);
    x0 = mod(x-2,edge_size)+1;
    x1 = mod(x,edge_size)+1;
    y0 = mod(y-2,edge_size)+1;
    y1 = mod(y,edge_size)+1;
    energy_j_old = -J*lattice(x,y)*(lattice(x0,y)+lattice(x1,y)+lattice(x,y0)+lattice(x,y1));
    energy_k_old = -K*lattice(x,y)*(lattice(x0,y0)+lattice(x0,y1)+lattice(x1,y1)+lattice(x1,y0));
    energy_j_new = -J*(-lattice(x,y))*(lattice(x0,y)+lattice(x1,y)+lattice(x,y0)+lattice(x,y1));
    energy_k_new = -J*(-lattice(x,y))*(lattice(x0,y)+lattice(x1,y)+lattice(x,y0)+lattice(x,y1));
    energy_new = energy_j_new + energy_k_new;
    energy_old = energy_j_old + energy_k_old;
    ratio = exp(-Beta*(energy_new - energy_old));
end

