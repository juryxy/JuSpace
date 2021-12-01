function [data_num] = cell2num_my(data)

        tt = isemptycell(data);
        data(tt==1)= {'NaN'};

% [r,c]=ind2sub(size(data), strmatch('NA', data, 'exact'));
m = size(data);
% 
% for l = 1:length(r)
% data{r(l),c(l)}='-4';
% end

data_num = zeros(m(1),m(2));

if isstr(data{1,1})
for i=1:m(1)
%     disp(i);
    for j=1:m(2)
        data_num(i,j) = str2double(data{i,j});
    end
end
else
  for i=1:m(1)
%     disp(i);
    for j=1:m(2)
        data_num(i,j) = data{i,j};
    end
  end
end