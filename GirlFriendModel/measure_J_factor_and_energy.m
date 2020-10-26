function [J_factor,energy] = measure_J_factor_and_energy(lattice,J,K)
     edge_size = length(lattice);
    energy = 0.0;
    J_factor = 0.0;
    for x_index = 1:1:edge_size
        for y_index = 1:1:edge_size
            x = x_index;
            y = y_index;
            x0 = mod(x-2,edge_size)+1;
            x1 = mod(x,edge_size)+1;
            y0 = mod(y-2,edge_size)+1;
            y1 = mod(y,edge_size)+1;
            energy_j = -J*lattice(x,y)*(lattice(x0,y)+lattice(x1,y)+lattice(x,y0)+lattice(x,y1))/2;
            energy_k = -K*lattice(x,y)*(lattice(x0,y0)+lattice(x0,y1)+lattice(x1,y1)+lattice(x1,y0))/2;
            J_factor = J_factor + lattice(x,y)*(lattice(x0,y)+lattice(x1,y)+lattice(x,y0)+lattice(x,y1))/2;
            energy = energy + energy_j + energy_k;
        end
    end 
    energy = energy / edge_size^2;
    J_factor = J_factor / edge_size^2;
end

