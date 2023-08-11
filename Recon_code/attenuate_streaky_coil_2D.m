function kSpace_attenuated = attenuate_streaky_coil_2D(kSpace,Angles)
[nx,np,nc]=size(squeeze(kSpace));
ray1 = -0.5:1/nx:0.5-1/nx;
for ii = 1:length(Angles)
    kdata2(:,ii) = ray1.*exp(1i*Angles(1,ii));
end
w = abs(real(kdata2) + 1i*imag(kdata2))./max(abs(real(kdata2(:))  + 1i*imag(kdata2(:)))); w(isnan(w))=1;
par.E = MCNUFFT_2D_GPU_kxyt_single( kdata2, w, ones(nx,nx,nc) );
par.y = squeeze(kSpace);
FBP_CoilSens = squeeze( par.E'*par.y );
FBP_CoilSens_norm = FBP_CoilSens./max(abs(FBP_CoilSens(:)));
%% manipulate signal intensity in k-space
[nx,ny,nc] = size( FBP_CoilSens_norm );
aa = [];      LocInfo = [];
for ii = 1:nc
    aa = abs(FBP_CoilSens_norm(:,:,ii));
    [MaxValue,idx] = max(aa(:));
    [x,y] = ind2sub(size(aa),idx);
    LocInfo(:,ii) = [x y  MaxValue];
end
clear aa
% % calculate the actual distances between the center of FOV and the max siganal intensity in each coil
SpatialRes = [1.5  1.5];  % spatial resolution in x, y, z, respectively
tmp = (LocInfo(1:end-1,:)-[nx/2+1;  ny/2+1]).*SpatialRes';
Distance = sqrt(sum( tmp.^2, 1) );
% % sort the coil distances from the center of FOV, based on the max signal intensity                
[val1,ind1] = sort( Distance );
% % divide into two coil groups: 1) close to heart; 2) far from heart
EndInd = round( nc/2 );
CloseCoilGroup = ind1(1:EndInd);
FarCoilGroup = ind1(EndInd+1:end);
% % find the max signal intensity within the close coil group
% % this will serve as the threshold for attenuation in the far coil group
MaxSignal_InCloseCoilGroup = max( LocInfo(end, CloseCoilGroup) ); 
% % sort the max signal intensities in all coils
[val2,ind2] = sort( LocInfo(end,:) );  
[val1'  ind1'  val2'  ind2'];
% % find a coil index which matchs to threshold in sorted data
% % the result must be the coils from the far coil group 
% % because the threshold is the max value in the close coil group. 
ind3 = find( val2 >= MaxSignal_InCloseCoilGroup );
% % determine the attenuating ratios for the selected coils in the far coil group 
% % which have higher signal intensities than the threshold
inds = ind2(ind3:end);
AttenuationRatio = MaxSignal_InCloseCoilGroup./val2(ind3:end);

% % attenuate
kSpace_attenuated = squeeze(kSpace);
for ii = 1:length(inds)
    sel_coil = inds(ii);
    kSpace_attenuated(:,:,sel_coil) = kSpace_attenuated(:,:,sel_coil).*AttenuationRatio(ii);
end
end