clc, clear all

%% add your path
path = [datapath patientname '/'];
files = dir( path );
%% make a list of the .dat files
count = 0;
for k = 3:size(files,1)
    tmp = lower( files(k).name );
    if strcmp( tmp(end-2:end), 'dat' )
        count = count + 1;
        filename_all{count} = [path files(k).name];
        disp( [num2str(count) ': ' files(k).name] );
    end
end
%% choose which file to reconstruct
filenum = 1;
disp( '*********************************************************' )
disp( ' Selected file name: ' )
disp( filename_all{filenum} );
disp( '*********************************************************' )
filename = filename_all{ filenum };
%% new method
dataraid = ReadSiemensMeasVD13_idea(filename);
i = size(dataraid,2);
% 1 for skyra,avanto, 2 for Aera
rawdata2 = dataraid(i).rawdata;
mdh = dataraid(i).loopcounters;
nMeasurements = length(dataraid(i).loopcounters);
for iM=1:nMeasurements
    cha(iM)        = mdh(iM,15);
    lin(iM)        = mdh(iM,1);
    slc(iM)        = mdh(iM,3);
    angle_mdb(iM)       = mdh(iM,39);
end
angle_mdb(1) = 0;
angle_mdb = angle_mdb./10000;
for iM=1:nMeasurements
     kSpace(:,lin(iM),slc(iM),cha(iM) ) = rawdata2(iM,:);   % rawdata2(:,iM) --> rawdata2(iM,:), vd13a with ReadSiemensMeasVD13_idea.m
     angles_all(lin(iM),slc(iM)) = angle_mdb(iM);
end

[nx,np,nz,nc] = size(kSpace);
clear rawdata2 cha lin slc angle_mdb
% tiny golden angle
N = 5; tau = (1.0+sqrt(5.0))/2.0;
golden_angle = pi / (tau+N-1);
Angles = mod((0:golden_angle:(np-1)*golden_angle),2*pi);
%% organize data
coil_index =1:nc;
%% parameter
nshare=24; 
nspokes= 9;
%%  iterative recon
for slice = 1:nz %1:nz 
display(['slice', num2str(slice)])

clear kdata1 kdata2
tic
kdata = kSpace(:,:,slice,coil_index);
kdata = attenuate_streaky_coil_2D(kdata,Angles);
kdata = coil_compression(kdata,8);
%% Low spatial resolution reconstruction. You can test for different matrix size
% clearvars -except basis kdata Angles slice nspokes img patientname recon_LRBW
MatrixSub=112; % 448 224 112 56

kdata1=kdata((nx-MatrixSub)/2+1:end-(nx-MatrixSub)/2,:,:); 
[xx,np,nc]=size(squeeze(kdata1));

kdata1 =  double(squeeze(kdata1));
ray1 = -0.5:1/xx:0.5-1/xx;
for ii = 1:length(Angles)
    kdata2(:,ii) = ray1.*exp(1i*Angles(1,ii));% position
end
w = abs(real(kdata2) + 1i*imag(kdata2))./max(abs(real(kdata2(:))  + 1i*imag(kdata2(:)))); w(isnan(w))=1;
par.E = MCNUFFT_2D_GPU_kxyt_single( kdata2, w, ones(xx,xx,nc) );
par.y = kdata1;
NUFFT_CoilSens = squeeze( par.E'*par.y );
[~,b1sens] = adapt_array_2d_st2( NUFFT_CoilSens );
b1sens = b1sens./max( abs(b1sens(:)) );
disp('b1 sensitivity done.')
par.E2 = MCNUFFT( kdata2, w, b1sens);
recon_timeaverage=par.E2'*par.y; %zero padded image %coil combined time averaged image
%% calculate mask +kwic
ny=nspokes+nshare;
nt=floor(np/nspokes);
ntt=floor(np/nspokes)-2;
%% no kwic and vs
%  w = ones(MatrixSub,ny,22);

%% subspace kwic
half_width=56;
WeightedContrast=ones(nx,ny,ntt);
ramp=(half_width)/(nshare/2);
for num=1:(nshare/2)-1
        WeightedContrast(169+ramp*(num-1):224,num,:)=eps;
        WeightedContrast(224:280-ramp*(num-1),num,:)=eps;
end
WeightedContrast(223:225,nshare/2,:)=eps; 
mm=fliplr(WeightedContrast);
WeightedContrast=WeightedContrast+mm-ones(nx,ny,ntt);
w=WeightedContrast(169:280,:,:); % subspace
% figure;
% imagesc(w(:,:,1));colormap(gray);
%% devide radial spokes
clear kdatau,clear ku, clear wu
for ii=1:nt-2
    kdatau(:,:,:,ii)=kdata1(:,(ii-1)*nspokes+1:(ii)*nspokes+nshare,:);
    ku(:,:,ii)= kdata2(:,(ii-1)*nspokes+1:(ii)*nspokes+nshare);
end
param.E = MCNUFFT_2D_GPU_kxyt(ku,w,b1sens);
param.y = kdatau;
recon_nufft=param.E'*param.y; %zero padded image
% figure;imagescn(abs((recon_nufft)),[],[],[],3);
%% CS in subspace
lambda = 0.02; % 
param.TV = TV_Temp();
param.nite = 9;
param.display = 1;
param.TVWeight = lambda*prctile(abs(recon_timeaverage(:)),95);
recon_cs = squeeze(recon_nufft);
record_time(slice,1) = toc;
ite = 3;
for n = 1 : ite
    recon_cs = CSL1NlCg(recon_cs,param);
end
record_time(slice,2) = toc;
% recon_LRes = (recon_cs)./max(abs(recon_cs(:)));
% recon_LRes = fliplr(recon_LRes);
% recon_LRes = flipud(recon_LRes);
% figure;imagescn(abs((recon_LRes(xx/4+1:xx/4*3,xx/4+1:xx/4*3,:))),[0 0.18],[],[],3);
%%
tmp=gather(abs(recon_cs));
[nx,ny,nt]=size(tmp);
Data_Seq=reshape(tmp,nx*ny,nt);
covariance=cov(Data_Seq);
[PC, V] = eig(covariance);
V = diag(V);
[junk, rindices] = sort(-1*V);
V = V(rindices);
basis = PC(:,rindices);
% figure,plot(V,'*')
%% load data
clear recon_cs param par recon_nufft kdata1 kdata2
K=4;% the number of basis, test for different values for your applications.  10

kdata1 = kdata;
kdata1 = double(squeeze(kdata1));
[nx,np,nc]=size(kdata1);
ray1 = -0.5:1/nx:0.5-1/nx;
for ii = 1:length(Angles)
    kdata2(:,ii) = ray1.*exp(1i*Angles(1,ii));
end
w = abs(real(kdata2) + 1i*imag(kdata2))./max(abs(real(kdata2(:))  + 1i*imag(kdata2(:)))); w(isnan(w))=1;
par.E = MCNUFFT_2D_GPU_kxyt_single( kdata2, w, ones(nx,nx,nc) );
par.y = kdata1;
NUFFT_CoilSens = squeeze( par.E'*par.y );
[dummy,b1sens] = adapt_array_2d_st2( NUFFT_CoilSens );
b1sens = b1sens./max( abs(b1sens(:)) );
% figure;imagescn(abs(fliplr(flipud(NUFFT_CoilSens(nx/4+1:nx/4*3,nx/4+1:nx/4*3,:,:)))),[0 0.00013],[],[],3);
% figure;imagescn(abs(fliplr(flipud(b1sens(nx/4+1:nx/4*3,nx/4+1:nx/4*3,:,:)))),[0 2],[],[],3);

disp('b1 sensitivity done.')
par.E2 = MCNUFFT( kdata2, w, b1sens);
recon_timeaverage=par.E2'*par.y; %zero padded image
%% calculate mask  kwic kspace
ny=nspokes+nshare;
nt=floor(np/nspokes);
%% no kwic
% w = ones(448,ny,22);
%%
% % % % %triangle kwic
half_width=56;
WeightedContrast=ones(nx,ny,ntt);
ramp=(half_width+3)/(nshare/2);
for num=1:(nshare/2)-1
        WeightedContrast(169+ramp*(num-1):224,num,:)=eps;
        WeightedContrast(224:280-ramp*(num-1),num,:)=eps;
        %WeightedContrast(256,nshare,:)=0;
end
WeightedContrast(223:225,nshare/2,:)=eps; 
%WeightedContrast(256,nshare/2,:)=eps;
mm=fliplr(WeightedContrast);
WeightedContrast=WeightedContrast+mm-ones(nx,ny,ntt);
w=WeightedContrast;
%% ring correction
kSpace_ring = IR_compensation(kdata1);
kdata_corr = ring_correction(kdata1,50:np,5);
%% CS full space
clear kdatau,clear ku, clear wu

for ii=1:nt-2
    kdatau(:,:,:,ii)=kdata1(:,(ii-1)*nspokes+1:(ii)*nspokes+nshare,:);
    ku(:,:,ii)= kdata_corr(:,(ii-1)*nspokes+1:(ii)*nspokes+nshare);
end
% param.E = MCNUFFT_2D_GPU_kxyt(ku,w,b1sens);
%kdatau = permute(kdatau,[1 2 4 3]);
param.y = kdatau;
% recon_cs = param.E'*param.y; %zero padded image
% figure;imagescn(abs(fliplr(flipud(recon_cs(nx/4+1:nx/4*3,nx/4+1:nx/4*3,:)))),[],[],[],3);
PCA=TempPCASub(basis(:,1:K));
param.PCA=PCA;
param.E=MCNUFFT_2D_GPU_kxyt_Sub(ku,w,b1sens,PCA);
recon_nufft=param.E'*param.y;
% figure;imagescn(abs(fliplr(flipud(recon_nufft(nx/4+1:nx/4*3,nx/4+1:nx/4*3,:)))),[],[],[],3);
%% CS
lambda = 0.02; % Hong, with PCA , with Pyo DCF
param.TV = TV_Temp();
param.nite = 9;
param.display = 1;
param.TVWeight = lambda*prctile(abs(recon_timeaverage(:)),95);
recon_cs = squeeze(recon_nufft);

ite = 3;
record_time(slice,3) = toc;
for n = 1 : ite
    recon_cs = CSL1NlCg_PCA_sub(recon_cs,param);
end
record_time(slice,4) = toc;
% figure;imagescn(abs(fliplr(flipud(recon_cs(nx/4+1:nx/4*3,nx/4+1:nx/4*3,:)))),[0 0.001],[],[],3);
tmp=abs(param.PCA'*(gather((single(gather(recon_cs/max(recon_cs(:))))))));
% figure;imagescn(abs(fliplr(flipud(tmp(nx/4+1:nx/4*3,nx/4+1:nx/4*3,:)))),[0 0.05],[],[],3);
tic
%% LRBW filter
param.block_size  = [8 8];
param.tau_size    = 0.2;
param.Iterations  = 2;
recon_LR = LRBW_Denoising_2Dt(tmp, param );
recon_time_2(slice,1)=toc;
recon_LRBW(:,:,:,slice) = recon_LR;%./prctile(abs(recon_LR),95);
recon(:,:,:,slice)=tmp;
end

%%
recon_pro=recon_LRBW;
recon_pro = (recon_pro)./max(abs(recon_pro(:)));
recon_pro = fliplr(recon_pro);
% recon_pro=flip(recon_pro,4);
recon_pro=flipud(recon_pro);


nx = size(recon_pro,1);
figure;imagescn(abs((recon_pro(nx/4+1:nx/4*3,nx/4+1:nx/4*3,:,7))),[0 0.2],[],[],3);
%%
%save your data

