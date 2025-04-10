function [files] = select_con_maps_forfMRI_my(work_dir,Folder_filt,con_map_filt,not_filt,frames)
% function [files] = select_con_maps_forfMRI_my(work_dir,Folder_filt,con_map_filt,not_filt,frames)
% not filter (not_filt) is optional
mg = cd;
cd(work_dir);

if exist('frames')
    files_all = cellstr(spm_select('ExtFPListRec',cd,con_map_filt,frames));
else
    files_all = cellstr(spm_select('FPListRec',cd,con_map_filt));
end
if exist('Folder_filt')==1
    if isempty(Folder_filt)
        Folder_filt = 'data';
    end
    if iscell(Folder_filt)
        for j = 1:length(Folder_filt)
            ind_all(:,j) = isemptycell(regexp(files_all,Folder_filt{j}))==0;            
        end
        files = files_all(sum(ind_all')'>0);
    else
    files = files_all(isemptycell(regexp(files_all,Folder_filt))==0);
    end
end

if exist('not_filt')==1
    if iscell(not_filt)
        for j = 1:length(not_filt)
            files = files(isemptycell(regexp(files,not_filt{j}))==1);
        end
    else
    files = files(isemptycell(regexp(files,not_filt))==1);
    end
end
cd(mg);

disp(char(files));


