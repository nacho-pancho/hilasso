function x_mask = mdlsMakeHann2D(wr)
    x_mask = zeros(wr,wr);
    wr2 = wr/2;
    h = hann(10*wr2);
    for i = 1:wr
        for j = 1:wr
            r = sqrt((i-0.5-wr2)^2+(j-wr2-0.5)^2);
            if r < wr2
                x_mask(i,j) = h(round(5*(wr2+r)));
            end
        end
    end
end
