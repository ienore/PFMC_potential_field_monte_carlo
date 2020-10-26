function energy = energy_radio_model(lattice,Kernel)
    Potential = conv2(lattice,Kernel,'same');
    energy = sum(lattice.*Potential,'all');
end

