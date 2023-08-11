% Li Feng, Ricardo Otazo, NYU, 2012
addpath('nufft_toolbox/');
clear all
addpath C:\Data\El\Hassan\functions\mapVBVD
addpath ('C:\Data\El\Hassan\functions\mapVBVD')
addpath ('C:\Data\El\Hassan\CS_Rawdata_Reconstruction')
addpath ('C:\Data\El\Hassan\Opening_raw_data\VD13\Files_for_VD13')
addpath ('C:\Data\El\Hassan\imagescn_R2008a')
addpath ('C:\Data\El\Hassan\mapVBVD')
addpath ('C:\Data\El\Hassan\BWLRcode')
addpath ('C:\Data\El\Hassan\Data_Northwestern')
path = ('C:\Data\El\Hassan\NorthWestern_studies\CAMRI-VIDIF\SaveData');
addpath ('C:\Data\El\Hassan')
%% nufft addpaths
addpath ('C:\Data\El\Hassan\irt\utilities')
addpath ('C:\Data\El\Hassan\irt\systems')
addpath ('C:\Data\El\Hassan\irt\nufft')
%%
%Using Siemens vb13 program to extract kpace data as an m file
filename = 'C:\Data\El\Hassan\NorthWestern_studies\2016_03_18_CAMRI-THMOM\Raw_data\meas_MID00116_FID19096_RT_Cine_Radial_Calib.dat';
dummy_on = 1;
twixdata = mapVBVD(filename);
twixdata{1,2}.image.flagRemoveOS = 0; %This line preserves oversampling in nx
kSpace=double(twixdata{1,2}.image());
timeframe = 1;
lambda1=0.05; %TTV weighting factor
%% VD13 load and organize here
timeframe = 1;
% kSpace=kdata;
kSpace = squeeze(kSpace);
    if(length(size(kSpace))==5)
        kSpace=permute(kSpace,[1 3 5 4 2]);
        [sr, sth, sz, sl, nc]=size(kSpace);
    elseif(length(size(kSpace))==4)
        kSpace=permute(kSpace,[1 3 4 2]);
        [sr, sth, sz, nc]=size(kSpace);
        sl=1;
        temp(:,:,:,1,:)=kSpace; clear kSpace
        kSpace = temp; clear temp
    end
%%  Main reconstruction code start here
[sr, sth, sz, sl, nc]=size(kSpace);
kSpace = kSpace(:,:,1:sz,:,:,:);
[sr,sth,sz,slc,lset] = size(double(kSpace));
 tic       
%% Main Slice Loop
addpath C:\Hassan\forGPU\Data_Northwestern
addpath C:\Data\El\Hassan\functions
addpath C:\Data\El\Hassan\BWLRcode
addpath C:\Data\El\Hassan\CS_Rawdata_Reconstruction
path = 'C:\Data\El\Hassan\NorthWestern_studies\2016_03_18_CAMRI-THMOM\SaveData\'
N=7;
clear recon_cs_TV;
         for slice_counter=1:slc
           display(num2str(slice_counter));
            Coil1=squeeze(double(kSpace(:,:,:,slice_counter,:)))*1e6;
            [mx,my,mz,mc,ms] = size(Coil1);
            Coil_temp = zeros(mx,my,mz,mc,ms);
            Coil_temp(:,:,1:4,:,:) = Coil1(:,:,mz-3:mz,:,:);
            Coil_temp(:,:,5:mz,:,:) = Coil1(:,:,1:mz-4,:,:);
            Coil1 = trajector_corr_v2_hemi(Coil_temp,2,5,pi/24,1);
%             keyboard
clear kdata2; clear mask2
for jj = 1:size(Coil1,4)
            [kdata2(:,:,:,jj),mask2(:,:,:,jj)]=rad_to_cart_bilinear_sparse_gold_fibo_mod(Coil1(:,:,:,jj),N);
end            
            if dummy_on == 1
                kdata2 = kdata2(:,:,35:sz-4,:);
                mask2 = mask2(:,:,35:sz-4,:);
            else
                kdata2 = kdata2(:,:,35:sz-4,:);
                mask2 = mask2(:,:,35:sz-4,:);
            end
            nx=sr; ny=sr; nt=sz; 
            %Coil maps generated here
            DC=squeeze(sum(kdata2,3));
            SM=squeeze(sum(mask2,3));
            SM(find(SM==0))=1;
            for ii=1:nc
               DC(:,:,ii)=DC(:,:,ii)./SM(:,:,ii);
            end
            %%% This ref is the unpaired FFT for the coil sensitivity maps. The
            %%% unpaired FFT for the kdata is withinEmat_2DPERF.mtimes
            ref=ifft2c_mri(DC); 
            [dummy,b1]=adapt_array_2d_st2(ref);
            b1=b1/max(b1(:));
         
            [nr,np,nt,nc] = size(Coil1);
            count = 1;
            for jj = 1:nt;
                for kk = 1:np
                    kdata_corrected(:,count,:) = squeeze(Coil1(:,kk,jj,:));
                    count = count +1;
                end
            end
ntviews = nt*np;
nspokes=12;
nx = size(kdata_corrected,1);
N = 7; tau = (1.0+sqrt(5.0))/2.0;
golden_angle = pi / (tau+N-1);
Angles = golden_angle*(0:nt*np-1);
Angles = mod(Angles,pi);
ray1 = -0.5:1/nx:0.5-1/nx;
for ii = 1:length(Angles)
    k(:,ii) = ray1.*exp(1i*Angles(ii));
end
% ramp filter 
w = abs(real(k) + 1i*imag(k))./max(abs(real(k(:))  + 1i*imag(k(:)))); w(isnan(w))=1;
% data dimensions
clear kdata;
for ch=1:nc,kdata_corrected(:,:,ch)=kdata_corrected(:,:,ch).*sqrt(w);end
nt=floor(ntviews/nspokes);
% crop the data according to the number of spokes per frame
k=k(:,1:nt*nspokes);
w=w(:,1:nt*nspokes);
% sort the data into a time-series
clear kdatau,clear ku, clear wu,clear kdatau_venc
for ii=1:nt
    kdatau(:,:,:,ii)=kdata_corrected(:,(ii-1)*nspokes+1:ii*nspokes,:);
    ku(:,:,ii)=k(:,(ii-1)*nspokes+1:ii*nspokes);
    wu(:,:,ii)=w(:,(ii-1)*nspokes+1:ii*nspokes);
end
ku = 1i.*real(ku)+imag(ku);
if dummy_on == 1
          ku_dummy = ku(:,:,35:size(ku,3));
          wu_dummy = wu(:,:,35:size(ku,3));
          kdatau_dummy = kdatau(:,:,:,35:size(ku,3));
else
          ku_dummy = ku(:,:,35:size(ku,3));
          wu_dummy = wu(:,:,35:size(ku,3));
          kdatau_dummy = kdatau(:,:,:,35:size(ku,3));
end
% multicoil NUFFT operator
param.E=MCNUFFT(ku_dummy,wu_dummy,b1);
% undersampled data
param.y=kdatau_dummy;
% clear kdata kdatau k ku wu w
% nufft recon
recon_nufft=param.E'*param.y;
% parameters for reconstruction
param.W = TV_Temp();
%% the L1 weight used
param.lambda=0.25*max(abs(recon_nufft(:))).*.45;
%%
param.nite = 10;
param.display=1;
fprintf('\n GRASP reconstruction \n')
tic
recon_cs=recon_nufft;
for n=1:4,
	recon_cs = CSL1NlCg_new(recon_cs,param);
%     save('('C:\Data\El\Hassan\NorthWestern_studies\2016_03_18_CAMRI-THMOM\Raw_data\save_data\recon_cs_lp025_25.mat','recon_cs')
end
recon_cs_TV(:,:,:,slice_counter) = recon_cs;
save('C:\Data\El\Hassan\NorthWestern_studies\2016_03_18_CAMRI-THMOM\SaveData\recon_cs_TV_so_far_2.mat','recon_cs_TV');
end
save('C:\Data\El\Hassan\NorthWestern_studies\2016_03_18_CAMRI-THMOM\SaveData\recon_cs_radial_calibration.mat','recon_csTV');



clear all;
path = 'C:\Data\El\Hassan\NorthWestern_studies\2016_03_18_CAMRI-THMOM\SaveData\'
file = 'C:\Data\El\Hassan\NorthWestern_studies\2016_03_18_CAMRI-THMOM\Raw_data\meas_MID00109_FID19089_RT_Cine_Cartesian_CS.dat';
filename = [file]
lambda1 = 0.01;  %TTV weighting factor
%% VB17 load and organize here
%% VD13 load and organize here
twixdata = mapVBVD2(filename);
% twixdata{1,2}.image.flagRemoveOS = 0; %This line preserves oversampling in nx
kSpace=twixdata{1,2}.image();
kSpace = squeeze(kSpace);
    if(length(size(kSpace))==5)
        kSpace=permute(kSpace,[1 3 5 4 2]); %nx ny nt sl nc
        [nx, ny, nt, sl, nc]=size(kSpace);
    elseif(length(size(kSpace))==4)
        kSpace=permute(kSpace,[1 3 4 2]);
        [nx, ny, nt, nc]=size(kSpace);
        sl=1;
    end
%% Main reconstruction code start here
tic    
kdata = kSpace;
[nx, ny, nt, sl, nc]=size(kdata);
kdata=permute(kdata,[1 2 3 5 4]); %Rearrange back from PCA
%% Main Slice Loop
for z = 1:sl
disp(num2str(z));
temp = kdata(:,:,:,:,z);
kdata2 = temp;
mask=zeros(size(kdata2,2),size(kdata2,3));
for zz=1:size(mask,2)
    mask(find(kdata2(size(kdata2,1)/2+1,:,zz,1)),zz)=1;
end
acc=size(kdata2,2)*size(kdata2,3)/sum(mask(:));
%raw data normalization
aux=ifft2c_mri(kdata2);sf=max(aux(:));aux=aux/sf;
kdata2=fft2c_mri(aux);
mask2 = zeros(nx,ny,nt);
for i = 1:nx
    mask2(i,:,:) = mask;
end
DC=squeeze(sum(kdata2,3));
SM=squeeze(sum(mask2,3));
SM(find(SM==0))=1;
for ii=1:nc
    DC(:,:,ii)=DC(:,:,ii)./SM;
end
%%% This ref is the unpaired FFT for the coil sensitivity maps. The
%%% unpaired FFT for the kdata is withinEmat_2DPERF.mtimes
ref=ifft2c_mri(DC); 
[dummy,b1]=adapt_array_2d_st2(ref);
b1=b1/max(b1(:));
param.W=TempFFT(3);param.L1Weight=lambda1/10;
param.TV = TV_Temp();param.TVWeight=lambda1;
param.nite = 9;
param.display = 1;
% iterative reconstruction ************************************************
param.E=Emat_2DPERF(mask2,b1); 
%recon the first image
param.y = kdata2;        
recon_zf=param.E'*param.y;                 %zeros filling reconstruction
recon_cs=recon_zf;
ite=4;
for n=1:ite,
	recon_cs = CSL1NlCg_Temp_TV(recon_cs,param);
end
recon_cs_TV(:,:,:,z)=recon_cs/max(recon_cs(:));

% not allowing reconstructions to be waisted, and furthermore
%allowing for reconstructions to pick up in case of failure
recon_cart_TV = recon_cs_TV;
filename = [path,'Cartesian_TTV_so_far.mat'];
save(filename,'recon_cart_TV');
end
toc
filename = [path,'Cartesian_TTV_recon_nx144_ny108_S12_rect_FOV_l01.mat'];
save(filename,'recon_cart_TV');












