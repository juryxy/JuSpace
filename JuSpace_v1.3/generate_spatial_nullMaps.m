function [data_permuted] = generate_spatial_nullMaps(atlas,data,N,opt_symmetry)
% [data_rand_weighted] = generate_spatial_nullMaps(atlas,data,N,opt_symmetry))
% atlas: file path to the atlas
% data: extracted data from JuSpace
% N: number of nullMaps
% opt_symmetry: 1 (default) generate symmetric null maps, else 0 -->
% asymetric

Y = spm_read_vols(spm_vol(atlas));
data_PET = data;

[a,b,c] = unique(Y(:));


if ~exist('opt_symmetry','var')
    opt_symmetry = 1;
end
if opt_symmetry == 1
    Y_half = zeros(size(Y));
    Y_half(1:floor(size(Y,1)./2),:,:) = Y(1:floor(size(Y,1)./2),:,:);
    [a1,b1,c1] = unique(Y_half(:));
end

% -----------------------
% compute center of mass per atlas region
for i = 2:length(a)
    Y_i = zeros(size(Y));
    Y_i(c==i) = 1;
    coord(i-1,:) = center_of_mass_my(Y_i);
end
% -----------------------
a1_ind_filt = find(a1~=0 & ~isnan(a1));
a1_filt = a1(a1_ind_filt);
if opt_symmetry == 1
    a_ind_filt = find(a~=0 & ~isnan(a));
    a_filt = a(a_ind_filt);
end
if opt_symmetry == 1
    for i = 1:length(a1_filt)
        ind_i = find(a1_filt(i)==a(2:end));
        coord_half(i,:) = coord(ind_i,:);
        data_fin(1,i) = data(ind_i);
    end
    coord_fin = coord_half;
else 
    data_fin = data;
    coord_fin = coord;
end

dist = [];
dist_vals = [];


for i = 1:length(coord_fin)
    for j = 1:length(coord_fin)
        dist(i,j) = sqrt(sum((coord_fin(i,:) - coord_fin(j,:)).^2)); % Euclidean distance for each pair
        dist_vals(i,j) = abs(data_fin(i) - data_fin(j)); % value distance for each pair
    end
end

[r] = corr(dist(:),dist_vals(:));
if r>0
    try
        parfor nn = 1:N
            disp(nn);
            data_rand_weighted = zeros(length(coord_fin),1);
            dist_valsr = zeros(length(coord_fin),length(coord_fin));
%             tic;
            rr = randperm(length(data_fin));
            data_rand = data_fin(rr);
            min_data = min(data_rand);
            max_data = max(data_rand);
                for i = 1:length(coord_fin)
                       for j = 1:length(coord_fin)
                                dist_valsr(i,j) = abs(data_rand(i) - data_rand(j));
                       end
                end
               [rr] = corr(dist(:),dist_valsr(:));
    %             rr = -1;
                std_dist = 1;
                % smooth data to induce spatial autocorrelation
                if rr<r
                    rr_prev = 0;
                    while rr<r && round(rr,4)>rr_prev
                        std_dist = std_dist+1;
                        for i = 1:length(coord_fin)
                            dist_i = dist(i,:);
                            dist_max = max(dist_i);
                            weight_i = exp(-dist_i.^2./(2.*std_dist)).*1./sqrt(2.*pi.*std_dist);
                            weight_i = weight_i./sum(weight_i);
                            data_rand_weighted(i,1) = data_rand*weight_i';
                        end

                        for i = 1:length(coord_fin)
                            for j = 1:length(coord_fin)
                                dist_valsr(i,j) = abs(data_rand_weighted(i) - data_rand_weighted(j));
                            end
                        end
                        rr_prev = round(rr,4);
                        [rr] = corr(dist(:),dist_valsr(:));
                    end
                else
                    data_rand_weighted = data_rand;
                end
            %rescale to original min and max
            min_nn = min(data_rand_weighted);
            data_rand_weighted_c = data_rand_weighted-min_nn;
            data_rand_weighted_c = (data_rand_weighted_c.*(max_data-min_data)./max(data_rand_weighted_c))+min_data;
            % rescale end

            %----------------------
            % Create a symmetric brain and project values to both hemispheres
            if opt_symmetry == 1
                Y_half_nn = zeros(size(Y));
                for ii = 1:length(a1_filt)
                    Y_half_nn(c1==ii+1) = data_rand_weighted_c(ii);
                end
                    Y2 = flip(Y_half_nn);
                    YY = cat(1,Y_half_nn(1:floor(size(Y,1)./2),:,:),Y2(floor(size(Y,1)./2)+1:end,:,:));

                for ii = 1:length(a_filt)
                    ind_ii = Y(:)==a_filt(ii);
                    ii_sel = YY(ind_ii);
                    data_rand_weighted_fin(ii,1) = mode(ii_sel(ii_sel~=0));
                end
             %----------------------
            else
                    data_rand_weighted_fin = data_rand_weighted_c;
            end
            data_permuted(nn,:) = data_rand_weighted_fin;
%             toc;
        end
    catch %if no parallel computing toolbox available
         for nn = 1:N
            disp(nn);
            data_rand_weighted = zeros(length(coord_fin),1);
            dist_valsr = zeros(length(coord_fin),length(coord_fin));
            tic;
            rr = randperm(length(data_fin));
            data_rand = data_fin(rr);
            min_data = min(data_rand);
            max_data = max(data_rand);
                for i = 1:length(coord_fin)
                       for j = 1:length(coord_fin)
                                dist_valsr(i,j) = abs(data_rand(i) - data_rand(j));
                       end
                end
               [rr] = corr(dist(:),dist_valsr(:));
              
    %             rr = -1;
                std_dist = 1;
                % smooth data to induce spatial autocorrelation
                if rr<r
                    rr_prev = 0;
                    while rr<r && round(rr,4)>rr_prev
                        std_dist = std_dist+1;
                        for i = 1:length(coord_fin)
                            dist_i = dist(i,:);
%                             dist_max = max(dist_i);
                            weight_i = exp(-dist_i.^2./(2.*std_dist)).*1./sqrt(2.*pi.*std_dist);
                            weight_i = weight_i./sum(weight_i);
                            data_rand_weighted(i,1) = data_rand*weight_i';
                        end

                        for i = 1:length(coord_fin)
                            for j = 1:length(coord_fin)
                                dist_valsr(i,j) = abs(data_rand_weighted(i) - data_rand_weighted(j));
                            end
                        end
                        rr_prev = round(rr,4);
                        [rr] = corr(dist(:),dist_valsr(:));
                         disp(rr);
                    end
                else
                    data_rand_weighted = data_rand;
                end
            %rescale to original min and max
            min_nn = min(data_rand_weighted);
            data_rand_weighted_c = data_rand_weighted-min_nn;
            data_rand_weighted_c = (data_rand_weighted_c.*(max_data-min_data)./max(data_rand_weighted_c))+min_data;
            % rescale end

            %----------------------
            % Create a symmetric brain and project values to both hemispheres
            if opt_symmetry == 1
                Y_half_nn = zeros(size(Y));
                for ii = 1:length(a1_filt)
                    Y_half_nn(c1==ii+1) = data_rand_weighted_c(ii);
                end
                    Y2 = flip(Y_half_nn);
                    YY = cat(1,Y_half_nn(1:floor(size(Y,1)./2),:,:),Y2(floor(size(Y,1)./2)+1:end,:,:));

                for ii = 1:length(a_filt)
                    ind_ii = Y(:)==a_filt(ii);
                    ii_sel = YY(ind_ii);
                    data_rand_weighted_fin(ii,1) = mode(ii_sel(ii_sel~=0));
                end
             %----------------------
            else
                    data_rand_weighted_fin = data_rand_weighted_c;
            end
            data_permuted(nn,:) = data_rand_weighted_fin;
            toc;
        end
    end
else
      
      disp('Spatial autocorrelation is zero or negative, no adjustment performed')
      for nn = 1:N
          rr = randperm(length(data));
          data_permuted(nn,:) = data(rr);
      end
end
