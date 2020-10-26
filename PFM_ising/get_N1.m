function N1 = get_N1(lattice,x,y)
    edge_size = length(lattice);
    x0 = mod(x-2,edge_size)+1;
    x1 = mod(x,edge_size)+1;
    y0 = mod(y-2,edge_size)+1;
    y1 = mod(y,edge_size)+1;
    N1 = lattice(x0,y)+lattice(x1,y)+lattice(x,y0)+lattice(x,y1);
end

