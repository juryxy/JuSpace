function [img3d] = resize_img_useTemp_imcalc(file,temp)
%-----------------------------------------------------------------------
% Job saved on 15-Jul-2015 09:33:17 by cfg_util (rev $Rev: 6134 $)
% spm SPM - SPM12 (6225)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------


[path] = fileparts(file);

file_temp = char(append_prefix_to_fileNames_my(temp,'rmyTemp'));
matlabbatch{1}.spm.util.imcalc.input = {
                                        temp
                                        file
                                        };
matlabbatch{1}.spm.util.imcalc.output = file_temp;
matlabbatch{1}.spm.util.imcalc.outdir = {path};
matlabbatch{1}.spm.util.imcalc.expression = 'i2';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 0;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
spm_jobman('run',matlabbatch);

img3d = spm_read_vols(spm_vol(file_temp));
file_temp = regexprep(file_temp,',1','');
delete(file_temp);
