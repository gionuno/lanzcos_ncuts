function v = lanczos(I,sI,sX,R,K,T)
    
    D = zeros(size(I,1),size(I,2));
    W = zeros(size(R,1),size(I,1),size(I,2));
    for x = 1:size(I,1)
    for y = 1:size(I,2)
        for r = 1:size(R,1)
            xr = x+R(r,1);
            yr = y+R(r,2);
            if 1 <= xr && xr <= size(I,1) && 1 <= yr && yr <= size(I,2)
                W(r,x,y) = exp(-0.5*sum((I(x,y,:)-I(xr,yr,:)).^2)/sI^2)*exp(-0.5*sum([R(r,1),R(r,2)].^2)/sX^2);
            end
        end
        D(x,y) = sum(W(:,x,y));
        W(:,x,y) = W(:,x,y) ;%/ D(x,y);
    end
    end
    
    v = randn(K,size(I,1),size(I,2));
    v(1,:,:) = v(1,:,:) / sqrt(sum(sum(v(1,:,:).^2)));
    DWv = randn(K,size(I,1),size(I,2));
    
    for t = 1:T
        disp(t);
        for x = 1:size(I,1)
        for y = 1:size(I,2)
            DWv(:,x,y) = sqrt(D(x,y))*v(:,x,y);
            for r = 1:size(R,1)
            	xr = x+R(r,1);
                yr = y+R(r,2);
                if 1 <= xr && xr <= size(I,1) && 1 <= yr && yr <= size(I,2)
                	DWv(:,x,y) = DWv(:,x,y) + W(r,x,y)*v(:,xr,yr) / sqrt(D(xr,yr));
                end
            end
            DWv(:,x,y) = DWv(:,x,y) / sqrt(D(x,y));
        end
        end
        H = zeros(K,K);
        for k = 2:K
            for l = 1:(k-1)
                H(k,l) = sum(sum(DWv(k,:,:) .* v(l,:,:)));
            end
        end
        for k = 2:K
            for l = 1:(k-1)
                DWv(k,:,:) = DWv(k,:,:) - H(k,l) * v(l,:,:);
            end
        end
        for k = 1:K
            v(k,:,:) = DWv(k,:,:) / sqrt(sum(sum(DWv(k,:,:).^2)));
        end
        
        %v = v / sqrt(sum(sum(v.^2)));
    end
end