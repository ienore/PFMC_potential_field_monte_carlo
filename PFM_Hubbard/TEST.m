NumInEdge  = 20;
site_index = 150;
m=5;
Slice = zeros([1,NumInEdge^2]);
[x,y] = IndexToCoor(site_index,NumInEdge);
sum_NN_m = 0.0;
for x_index = 1:1:NumInEdge
    delta_x = min([abs(x_index-x),NumInEdge-abs(x_index-x)]);
    for y_index = 1:1:NumInEdge
        delta_y = min([abs(y_index-y),NumInEdge-abs(y_index-y)]);
        if delta_x + delta_y == m
            index = CoorToIndex(x_index,y_index,NumInEdge);
            Slice(index) = 1;
        end
    end
end     
show = zeros([NumInEdge,NumInEdge]);
for x_index = 1:1:NumInEdge
    for y_index = 1:1:NumInEdge
        index = CoorToIndex(x_index,y_index,NumInEdge);
        show(x_index,y_index) = Slice(index);
    end
end
pcolor(show);