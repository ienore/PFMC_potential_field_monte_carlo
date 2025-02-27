function [J_factor,C_factor,energy] = measure_C_factor_and_energy(lattice,connection,J,K)
    max_friend = length(connection(1,:))/2;
    edge_size = length(lattice);
    energy = 0.0;
    C_factor = 0.0;
    J_factor = 0.0;
    for x_index = 1:1:edge_size
        for y_index = 1:1:edge_size
            x = x_index;
            y = y_index;
            x0 = mod(x-2,edge_size)+1;
            x1 = mod(x,edge_size)+1;
            y0 = mod(y-2,edge_size)+1;
            y1 = mod(y,edge_size)+1;
            energy_j = -J*lattice(x,y)*(lattice(x0,y)+lattice(x1,y)+lattice(x,y0)+lattice(x,y1));
            energy_k = -K*lattice(x,y)*(lattice(x0,y0)+lattice(x0,y1)+lattice(x1,y1)+lattice(x1,y0));
            connection_index = (x_index - 1)*edge_size + y_index;
            for friend_index = 1:1:max_friend
                if connection(connection_index,2*(friend_index-1)+1) ~= 0 && connection(connection_index,2*(friend_index-1)+2) ~= 0
                    aim_x = connection(connection_index,2*(friend_index-1)+1);
                    aim_y = connection(connection_index,2*(friend_index-1)+2);
                    C_factor = C_factor + lattice(x,y)*lattice(aim_x,aim_y);
                end
            end
            J_factor = J_factor + lattice(x,y)*(lattice(x0,y)+lattice(x1,y)+lattice(x,y0)+lattice(x,y1));
            energy = energy + energy_j + energy_k;
        end
    end 
    
end

