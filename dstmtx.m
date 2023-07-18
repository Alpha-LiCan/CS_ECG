function ww = dstmtx(N)
    ww = zeros(N);
    
    for i = 1:N
        for j = 1:N
            ww(i, j) = sqrt(2/N) * sin((i-0.5) * (j-0.5) * pi/N);
        end
    end
    
    ww(1, :) = ww(1, :) / sqrt(2);
end