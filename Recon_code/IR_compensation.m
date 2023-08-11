function kSpace_comp = IR_compensation(kSpace)
[nx,np,nc] = size(kSpace);
kSpace_intensity = squeeze(sum(kSpace(nx/2-1:nx/2+1,:,:),1));
k_int = mean(abs(kSpace_intensity(:,:)),2);
k_int = k_int./max(k_int(:));
[~,inv_point] = min(k_int(:));
x1 = (1:inv_point-1)'; 
% y1 = k_int(1:inv_point-1);
x2 = (inv_point:np)';
x22 = (1:np-inv_point+1)';
y2 = k_int(inv_point:np);
g = fittype('a-b*exp(-c*x)');
% f1 = fit(x1,y1,g,'StartPoint',[[ones(size(x1)), -exp(-x1)]\y1; 1]);
f2 = fit(x22,y2,g,'StartPoint',[[ones(size(x22)), -exp(-x22)]\y2; 1]);
% figure; 
% hold on
% plot(1:np,k_int,'b-','linewidth',2);
% plot(x1,f1(x1),'r--','linewidth',2);
% plot(x2,f2(x22),'y--','linewidth',2);
% hold off

% fit_int(1:length(x1)) = f1(x1);
fit_int = ones(1,np);
fit_int(length(x1)+1:np) = f2(x2);% k_int(length(x1)+1:np);%f2(x2);
for i = 1:nc
    kSpace_comp(:,:,i) = kSpace(:,:,i)./fit_int;
end
% figure;imagescn(abs((kSpace)),[],[],[],3);
% figure;imagescn(abs((kSpace_comp)),[],[],[],3);
end