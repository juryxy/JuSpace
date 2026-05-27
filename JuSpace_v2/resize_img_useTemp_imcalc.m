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
% Uses nearest-neighbor assignment.


file_check = contains(file,'.nii,');
temp_check = contains(template,'.nii,');

if file_check
file = strrep(file,',1','');
end

if temp_check
template = strrep(template,',1','');
end



    srcInfo = niftiinfo(file);
    tmplInfo = niftiinfo(template);

    srcImg = niftiread(srcInfo);

    srcMat  = local_nifti_affine(srcInfo);
    tmplMat = local_nifti_affine(tmplInfo);

    srcDim = size(srcImg);
    outDim = double(tmplInfo.ImageSize(1:3));

    if numel(srcDim) < 3
        srcDim(3) = 1;
    end
    if numel(outDim) < 3
        outDim(3) = 1;
    end

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

    % Map template voxel coordinates to source voxel coordinates.
    Pin = srcMat \ (tmplMat * Pout);

    Xi = reshape(Pin(1,:), outDim);
    Yi = reshape(Pin(2,:), outDim);
    Zi = reshape(Pin(3,:), outDim);

    % Nearest-neighbor voxel assignment
    Xi = round(Xi);
    Yi = round(Yi);
    Zi = round(Zi);

    img3d = zeros(outDim, 'like', srcImg);

    valid = ...
        Xi >= 1 & Xi <= srcDim(1) & ...
        Yi >= 1 & Yi <= srcDim(2) & ...
        Zi >= 1 & Zi <= srcDim(3);

    srcIdx = sub2ind( ...
        srcDim(1:3), ...
        Xi(valid), ...
        Yi(valid), ...
        Zi(valid));

    img3d(valid) = srcImg(srcIdx);
end


function M = local_nifti_affine(info)
%LOCAL_NIFTI_AFFINE Return 4x4 voxel-to-world affine from niftiinfo output.
%
% Returns affine in column-vector convention:
%
%   world = M * [i; j; k; 1]
%
% using MATLAB/SPM-style 1-based voxel coordinates.

    if isfield(info, 'Transform') && ~isempty(info.Transform)
        T = info.Transform.T';

        % Convert MATLAB/NIfTI 0-based voxel transform to 1-based voxel coords.
        shift = eye(4);
        shift(1:3,4) = -1;

        M = T * shift;
    else
        pixdim = double(info.PixelDimensions(1:3));

        M = eye(4);
        M(1,1) = pixdim(1);
        M(2,2) = pixdim(2);
        M(3,3) = pixdim(3);
    end
end