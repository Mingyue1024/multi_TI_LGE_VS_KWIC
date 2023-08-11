function kSpace_comp = coil_compression(kSpace,N)
if length(size(kSpace))==4
[nx,ny,nz,nc] = size(kSpace);
    if(nc>N)
        no_comp=N;
    else
        no_comp=nc;
    end
data = reshape(kSpace,[nx*ny*nz nc]);
[coeff,~,~] = pca(data);
compressed_data = data*coeff(:,1:no_comp);
kSpace_comp = reshape(compressed_data,[nx,ny,nz,no_comp]);
else %% length(size(kSpace))==3
[nx,ny,nc] = size(kSpace);
    if(nc>N)
        no_comp=N;
    else
        no_comp=nc;
    end
data = reshape(kSpace,[nx*ny nc]);
[coeff,~,~] = pca(data);
compressed_data = data*coeff(:,1:no_comp);
kSpace_comp = reshape(compressed_data,[nx,ny,no_comp]);    
end
end