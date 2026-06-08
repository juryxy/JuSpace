function [img3d] = resize_img_useTemp_imcalc(file,temp)

    if nargin < 2
        error('Usage: img3d = resize_img_useTemp_imcalc(file,temp)');
    end

    Vfile = spm_vol(file);
    Vtemp = spm_vol(temp);

    img3d = zeros(Vtemp.dim, 'double');

    for z = 1:Vtemp.dim(3)
        M = inv(Vfile.mat) * Vtemp.mat * spm_matrix([0 0 z]);
        img3d(:,:,z) = spm_slice_vol(Vfile, M, Vtemp.dim(1:2), 0);
    end

end