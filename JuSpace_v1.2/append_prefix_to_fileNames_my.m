function [files_new] = append_prefix_to_fileNames_my(files,prefix)
files = cellstr(files);

files_new = files;
for i = 1:length(files)

    [path,name,ext] = fileparts(files{i});

    name_new = [prefix name];
    files_new{i} = fullfile(path,[name_new ext]);
end
    
    