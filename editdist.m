%edit distance calculation
function edr=editdist(X,Y,e)
    [si_x,k] = size(X);
    [si_y,k] = size(Y);
    D = zeros(si_x,si_y);
    D(1,1) = 0;
    for i=2:si_x
        D(i,1) = i-1;
    end
    
    for i=2:si_y
        D(1,i) = i-1;
    end
    
    for i=2:si_x
       for j=2:si_y
            m = X(i,:);
            n = Y(j,:);
            val = 0;
            if(dist(m,n) > e)
                val = 1;
            end
            v = [D(i,j-1)+1 D(i-1,j)+1 val+D(i-1,j-1)];
            D(i,j) = min(v);
       end
    end
    edr = D(si_x,si_y);
end
