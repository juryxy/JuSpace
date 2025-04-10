function [img3d] = resize_img_useTemp_imcalc(file, template)
    % Load input and template NIfTI images

    inputVol = spm_read_vols(spm_vol(file));
    templateVol = spm_read_vols(spm_vol(template));

    % Determine the target size (dimensions) from template
    targetSize = size(templateVol);

    % Resize the input volume to match the template size
    [img3d]  = resize3Dvolume(inputVol, targetSize);
    disp('Resize successful')
end


function resizedVol = resize3Dvolume(vol, targetSize)
    % vol: 3D input array
    % targetSize: [nx, ny, nz] new size
    
    [X, Y, Z] = ndgrid(...
        linspace(1, size(vol,1), targetSize(1)), ...
        linspace(1, size(vol,2), targetSize(2)), ...
        linspace(1, size(vol,3), targetSize(3)));
    
    resizedVol = interp3(vol, Y, X, Z, 'nearest', 0);  % 'linear' or 'nearest'
end