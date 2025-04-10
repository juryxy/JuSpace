function coord = center_of_mass_my(Y_b)

m = size(Y_b);

l = zeros(m(1),m(2));
n = zeros(m(1),m(3));
gg = round(m(3)./4);
for i = 1:m(3)
    l = l+Y_b(:,:,i);
end

[maxint,ind_max_y] = max(mean(l));
[maxint,ind_max_x] = max(mean(l'));


for i = 1:m(1)
    dd(:,:) = Y_b(:,i,:);
    n = n+dd;
end
% [a,b] = max(n(:));
% [ind_max_z,b]=ind2sub([m(1),m(3)],b);
% % 
[Y,I] = max(n);
[lg,ind_max_z] = max(Y(gg:end));
ind_max_z = ind_max_z + gg;
coord = [ind_max_x, ind_max_y,ind_max_z];