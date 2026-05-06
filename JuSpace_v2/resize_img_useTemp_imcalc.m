function [img3d]  = resize_img_useTemp_imcalc(file, template)
%SPM_LIKE_RESIZE_IMG Resize/resample NIfTI file onto template NIfTI grid.
%
% Inputs
%   file     : path to source .nii file
%   template : path to template .nii file
%
% Output
%   img3d    : source image resampled onto template voxel grid
%
% Requires MATLAB's niftiinfo/niftiread. Uses no SPM internal functions.

    % Read source and template headers/images
    srcInfo = niftiinfo(file);
    tmplInfo = niftiinfo(template);

    srcImg = double(niftiread(srcInfo));

    % Build voxel-to-world affine matrices
    srcMat  = local_nifti_affine(srcInfo);
    tmplMat = local_nifti_affine(tmplInfo);

    srcDim = size(srcImg);
    outDim = double(tmplInfo.ImageSize(1:3));

    % Handle possible 2D NIfTI edge case
    if numel(srcDim) < 3
        srcDim(3) = 1;
    end
    if numel(outDim) < 3
        outDim(3) = 1;
    end

    % Output voxel grid in MATLAB/SPM-style 1-based voxel coordinates
    [Xo, Yo, Zo] = ndgrid( ...
        1:outDim(1), ...
        1:outDim(2), ...
        1:outDim(3));

    Pout = [
        Xo(:)';
        Yo(:)';
        Zo(:)';
        ones(1, numel(Xo))
    ];

    % Core SPM-equivalent mapping:
    %
    % input_voxel = inv(srcMat) * tmplMat * output_voxel
    %
    % This maps every template voxel to the corresponding floating-point
    % source voxel coordinate, including origin, orientation, voxel size,
    % translation, and flips encoded in the NIfTI affine.
    Pin = srcMat \ (tmplMat * Pout);

    Xi = reshape(Pin(1,:), outDim);
    Yi = reshape(Pin(2,:), outDim);
    Zi = reshape(Pin(3,:), outDim);

    % Input grid
    [Xin, Yin, Zin] = ndgrid( ...
        1:srcDim(1), ...
        1:srcDim(2), ...
        1:srcDim(3));

    % Interpolate source image onto template grid
    img3d = interpn( ...
        Xin, Yin, Zin, ...
        srcImg, ...
        Xi, Yi, Zi, ...
        'linear', ...
        0);
end


function M = local_nifti_affine(info)
%LOCAL_NIFTI_AFFINE Return 4x4 voxel-to-world affine from niftiinfo output.
%
% MATLAB stores NIfTI transforms as affine3d objects where the matrix is
% commonly arranged for row-vector multiplication. This helper converts it
% to the column-vector convention:
%
%   world = M * [i; j; k; 1]
%
% with 1-based voxel coordinates, matching the SPM-style mapping above.

    if isfield(info, 'Transform') && ~isempty(info.Transform)
        T = info.Transform.T';

        % MATLAB's nifti transform is usually 0-based voxel-index oriented.
        % Convert to 1-based voxel coordinates:
        %
        % world = T * ([i; j; k; 1] - [1; 1; 1; 0])
        %
        shift = eye(4);
        shift(1:3,4) = -1;

        M = T * shift;

    else
        % Fallback: construct a simple affine from PixelDimensions.
        % This is less informative because it has no true scanner-space
        % origin or rotation.
        pixdim = double(info.PixelDimensions(1:3));

        M = eye(4);
        M(1,1) = pixdim(1);
        M(2,2) = pixdim(2);
        M(3,3) = pixdim(3);
    end
end