I = double(imread('lena.jpg'))/255.0;
I = I(1:4:end,1:4:end,:);

EI  = mean(mean(mean(I)));
EI2 = mean(mean(mean((I-EI).^2))); 

sI = 1.0*sqrt(EI2);
sX = 5.0;

r = 4;
R = [];
for a = (-r):r
for b = (-r):r
    if a ~= 0 || b ~= 0
        R = [R;[a,b]];
    end
end
end


T = 100;
K = 5;

v = lanczos(I,sI,sX,R,K*K,T);

figure;
for k = 1:K
    for l = 1:K
        subplot(K,K,K*(k-1)+l);
        aux = reshape(v(K*(k-1)+l,:,:),[size(I,1),size(I,2)]);
        aux = aux - min(min(aux));
        aux = aux / max(max(aux));
        imshow(aux);
        colormap jet;
    end
    
end

