function lrbw = LRBW_Denoising_2Dt(input_img, param )
% % Example)
% % param.block_size  = [4 4];
% % param.tau_size    = 0.04;
% % param.Iterations  = 1;
% % lrbw = LRBW_Denoising_2Dt( input_img, param )
% % % % % tau = tau_size*max( abs(recon_cs(:)) ); % set a threshold

fprintf( 'Low-rank Block-wise denoising for temporal 2D...\n' )


bloc_x = param.block_size(1);
bloc_y = param.block_size(2);
[NX NY NT] = size( input_img );
tau = param.tau_size * max( abs(input_img(:)) ); % Threshold 
temp_new_ga = padarray(input_img,[bloc_x-1  bloc_y-1],'symmetric','both');

lr_term_update = zeros( size(input_img) );
for iter = 1:param.Iterations
% %     fprintf( 'Iterations:  %d  out of  %d\n', iter,  param.Iterations )
    parfor iX = 1:NX  %/bloc_x
        for iY = 1:NY  %/bloc_y
            Bx = iX + bloc_x - 1;
            By = iY + bloc_y - 1;
% %             bloc = squeeze( temp_new_ga(Bx-(bloc_x-1):Bx+(bloc_x-1),By-(bloc_y-1):By+(bloc_y-1), :) );
            bloc = squeeze( temp_new_ga( Bx-(bloc_x-1):Bx+(bloc_x-1), By-(bloc_y-1):By+(bloc_y-1), :) );
            tempggg = reshape( bloc,[ size(bloc,1)*size(bloc,2)   size(bloc,3)] );
            if size(bloc,1)*size(bloc,2)  <  size(bloc,3)
                tempggg = tempggg';
            end
            
            [U,S,V] = svd( tempggg, 0 );
            S = diag(S);
            S = (S-tau);
            tmp_ind = find( S < 0 ); 
            S(tmp_ind) = 0;
            tempggg = U*diag(S)*V';
            if size(bloc,1)*size(bloc,2)  <  size(bloc,3)
                tempggg = tempggg';
            end
            tempggg = reshape( tempggg, [size(bloc,1)   size(bloc,2)  size(bloc,3)] );            
            lr_term_update( iX, iY, : ) = tempggg( 1+(bloc_x-1),1+(bloc_y-1), : );
        end
    end
end
lrbw = lr_term_update;
