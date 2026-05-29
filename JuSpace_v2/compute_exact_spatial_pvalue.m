function [p_exact,dist_rand,stats] = compute_exact_spatial_pvalue(D1,data_PET,atlas,res,Nperm,options,filesPET,T1,stats,d)
% function [p_exact,dist_rand,stats] = compute_exact_spatial_pvalue(D1,D2,data_PET,res,Nperm,options,filesPET,T1,stats,d)
% the functions allow to compute exact p-values based on the output from
% the compute_DomainGauges function, see help for this function for the
% input parameters:
% D1: Data for spatial correlation
% data_PET: extracted PET data
% atlas: full path to the atlas file
% res: Results for spatial correlation as provided by Results.res in
% compute_DomainGauges
% Nperm: Number of permutations
% options: options from the JuSpace toolbox
% filesPET: cell array with full path to PET files
% link to the T1 template if selected
% stats: stats results as provided by Results.stats
% d: optional messenger handle to update dialog window
% options: a numeric array, i.e. [1 1]
% first index indicates the computing option
    % option(1) = 1 --> es between
    % option(1) = 2 --> es within
    % option(1) = 3 --> mean list 1
    % option(1) = 4 --> list 1 each
    % option(1) = 5 --> ind z-score list 1 to list 2
    % option(1) = 6 --> pair-wise difference list 1 to list 2
    % option(1) = 7 --> leave one out from list 1
    % option(1) = 8 --> list 1 each compares against null distribution of
    % correlation coefficients
% second index indicates the analysis option
    % option(2) = 1 --> % Spearman correlation
    % option(2) = 2 --> % Pearson correlation
    % option(2) = 3 --> % multiple linear regresion
% IMPORTANT: exact spatial p-value is only supported for options(1) = 3, 4
% and 8

% if Nperm>1000
% catch

if isdeployed
     [~, ~] = system('path');
     path_ju = pwd; 
else
     path_ju = fileparts(which('JuSpace'));
end

[dd,atlas_name] = fileparts(atlas);
path_maps = fullfile(path_ju,'nullMaps',atlas_name);

if not(isfolder(path_maps))
    mkdir(path_maps)
end

 for i = 1:length(filesPET)
    [dd,PET_name] = fileparts(filesPET{i});
    null_path = fullfile(path_maps,[PET_name '.mat']);
    % opts.outputDirectory = path_maps;
    % opts.dataMapName = filesPET{i};
    % opts.useParallel = 1;
    if exist(null_path,'file')
        data_permuted = load(null_path);
        if size(data_permuted.data_permuted,1)< Nperm
            [data_permuted_n] = generate_spatial_nullMaps(atlas,data_PET(i,:),Nperm-size(data_permuted.data_permuted,1),1);
            data_permuted = [data_permuted.data_permuted; data_permuted_n];
            data_perm_all{i} = data_permuted;
            save(null_path,'data_permuted');
        else
            data_perm_all{i} = data_permuted.data_permuted;
        end
    else
%         for i = 1:size(data_PET,1)
        disp(['Generating permuted maps for PET map for ' PET_name]);
        if exist('d','var')
            d.Message = ['Generating permuted maps for PET map for ' PET_name];
        end
        [data_permuted] = generate_spatial_nullMaps(atlas,data_PET(i,:),Nperm,1);
        data_perm_all{i} = data_permuted;
        save(null_path,'data_permuted');
%         end
        disp('Generating PET maps completed');
    end
end

switch options(1)
    case 3
        if size(D1,1)==1
            data = D1;
        else
            data = mean(D1);
        end
    case {4,8}
         data = D1;
    case 7
         for nn = 1:size(D1,1)
             data(nn,:) = (D1(nn,:)-mean(D1([1:nn-1 nn+1:end],:)))./std(D1([1:nn-1 nn+1:end],:));
         end
end


    

switch options(2)

    case 2
        for i = 1:size(data_PET,1)

            if options(4) == 1
                data_ij = removenan_my([data',data_perm_all{i}',T1']);
                r_i = partialcorr(data_ij(:,1:size(data',2)),data_ij(:,size(data',2)+1:end-1),data_ij(:,end));
            else
                data_ij = removenan_my([data',data_perm_all{i}']);
                r_i = corr(data_ij(:,1:size(data',2)),data_ij(:,size(data',2)+1:end));
            end

            if options(1) == 8
                r_all{i} = mean(r_i);
            else
                r_all{i} = r_i;
            end
        end
    case 1
        for i = 1:size(data_PET,1)
%             r_i = corr(data',data_perm_all{i}','type','Spearman');
            if options(4) == 1
                data_ij = removenan_my([data',data_perm_all{i}',T1']);
                r_i = partialcorr(data_ij(:,1:size(data',2)),data_ij(:,size(data',2)+1:end-1),data_ij(:,end),'type','Spearman');
            else
                data_ij = removenan_my([data',data_perm_all{i}']);
                r_i = corr(data_ij(:,1:size(data',2)),data_ij(:,size(data',2)+1:end),'type','Spearman');
            end

            if options(1) == 8
                r_all{i} = mean(r_i);
            else
                r_all{i} = r_i;
            end
        end
    case 3

        for j = 1:Nperm
            x = data';
            Y = [];
            for k = 1:size(data_PET,1)
                Y = [Y data_perm_all{k}(j,:)'];
            end
            if options(4)==1
                Y = [Y T1'];
            end
            data_ij = removenan_my([x Y]);

            x = data_ij(:,1:size(data',2));
            x = zscore(x); %z-score
            
            Y = data_ij(:,size(data',2)+1:end);
            Y = zscore(Y); 
            
            Y = [ones(length(x),1) Y];
            r_i = [Y\x]';
            
            for k = 1:size(x,2)
                yhat = Y*r_i(k,:)';
                ybar = mean(x(:,k));
                residuals = x(:,k) - yhat;
                sse = norm(residuals)^2;    % sum of squared errors
                sst = norm(x(:,k) - ybar)^2;     % total sum of squares;
                rsquare_all(j,k) = 1 - sse ./ sst;    
            end
            
            if options(4)==1
                r_i = r_i(:,2:end-1);
            else
                r_i = r_i(:,2:end);
            end
%                 if options(1) == 8
%                     r_all_j(j,:) = mean(r_i);
%                 else
                r_all_j{j,1} = r_i;
                
         
            
%                 end
        end

        if options(1) == 8
            for k = 1:size(data_PET,1)
                for j= 1:Nperm
                    r_all_j_k{k}(:,j) = r_all_j{j}(:,k);
                end
            end
            for k = 1:size(data_PET,1)
                r_all{k} = mean(r_all_j_k{k});
            end
            
        else
            for k = 1:size(data_PET,1)
                for j= 1:Nperm
                    r_all{k}(:,j) = r_all_j{j}(:,k);
                end
            end
        end

            
end
    



    for i = 1:size(res,1)
        
        if options(2) == 3
            p_exact_f(i,1) = 1-(sum(stats.rsquare(i)>=rsquare_all(:,i))+1)./(length(rsquare_all)+1);
        else
            p_exact_f = '';
        end
        
        for j = 1:size(res,2)
            r_i_j = r_all{j}(i,:)';
            if isnan(res(i,j)) | isnan(r_i_j)
                p_exact(i,j) = NaN;
            else
                p_exact(i,j) = (sum(abs(r_i_j)>=abs(res(i,j)))+1)./(length(r_i_j)+1);
            end
        end
    end
    stats.p_exact_f = p_exact_f;
    % p_exact = (sum(abs(r_all)>=abs(res))+1)./(length(r_all)+1);
    p_exact = p_exact';
    p_exact = p_exact(:)';
    dist_rand = r_all;
