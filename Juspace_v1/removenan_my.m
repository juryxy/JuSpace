function [filt_dat,ind_filt] = removenan_my(data,col_ids)
% function [filt_dat,ind_filt] = removenan_my(data,col_ids)
if ~exist('col_ids','var')
    col_ids = ':';
end


if iscell(data)
m = [];
data_red = data(:,col_ids);
for i = 1:size(data_red,1)
    for j = 1:size(data_red,2)
        if isequal(data_red(i,j),{'NaN'})==1
            m = [m; i];
        end
    end
end
ff = unique(m);   
aa = [1:size(data_red,1)]';

gg = setdiff(aa,ff);
filt_dat = data(gg,:);

else
    data_red = data(:,col_ids);
    ind_filt = find(sum(isnan(data_red'),1)==0)';
    filt_dat = data(ind_filt,:);
end