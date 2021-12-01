function [data_num] = num2cell_my(data)

m = size(data);
try 
    data = cell2num_my(data);
end

for i=1:m(2)
        data_num(:,i) = cellstr(num2str(data(:,i)));
end

if m(2)<1
    data_num = cell(0,1);
end