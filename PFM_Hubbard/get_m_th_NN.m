function sum_NN_m = get_m_th_NN(Sigma,NumInEdge,site_index,time_index,m)
    Slice = Sigma(time_index,:);
    [x,y] = IndexToCoor(site_index,NumInEdge);
    sum_NN_m = 0.0;
    for x_index = 1:1:NumInEdge
        delta_x = min([abs(x_index-x),NumInEdge-abs(x_index-x)]);
        for y_index = 1:1:NumInEdge
            delta_y = min([abs(y_index-y),NumInEdge-abs(y_index-y)]);
            if delta_x + delta_y == m
                sum_NN_m = sum_NN_m + Slice(x_index,y_index);
            end
        end
    end       
end

