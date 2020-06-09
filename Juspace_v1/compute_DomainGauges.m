function [res,p_all,stats,data, D1,D2,data_PET,Resh,T1] = compute_DomainGauges(list1,list2,files_PET,atlas, options,image_save)
% [res,p_all,stats,data, D1,D2,data_PET] = compute_DomainGauges(list1,list2,files_PET,atlas, options, image_save)
% Inputs:
% list1, list2, files_PET are cellarrays of strings containing the filepath
% atlas cellstr of filepath to the atlas to use for DomainGauges, i.e.
% /atlas/m_labels_Neuromorphometrics.nii
% options: a numeric array, i.e. [1 1]
% first index indicates the computing option
    % option(1) = 1 --> es between
    % option(1) = 2 --> es within
    % option(1) = 3 --> mean list 1
    % option(1) = 4 --> list 1 each
    % option(1) = 5 --> ind z-score list 1 to list 2
    % option(1) = 6 --> pair-wise difference list 1 to list 2
    % options(1) = 7 --> leave one out from list 1
    % options(1) = 8 --> list 1 each compares against null distribution of
    % correlation coefficients
% second index indicates the analysis option
    % option(2) = 1 --> % Spearman correlation
    % option(2) = 2 --> % Pearson correlation
    % option(2) = 3 --> % multiple linear regresion
% if option(1) < 4 --> image_save: filepath (char) to save the produced image
% Outputs:
% [res,p_all,stats,data, D1,D2,data_PET]
% res --> Fisher's z transformed correlation / multiple linear regression coefficient matrix, rows correspond to
% input files, columns to PET maps
% p_all --> p-values from one-sample t-test testing if the Fisher's z
% transformed correlations /regression coefficients are different from zero
% stats --> Summary statistics
% stats.CorrOrig --> not Fisher's z transformed individual correlation
% coefficients
% stats.p_ind --> p-values for correlations (CorrOrig) of individual files with PET
% maps 
% stats.res_ind --> Individual Fisher's z-transformed correlation
% coefficients / regression coefficients
% stats.ci95 --> 95% confidence interval of the one sample t-test for
% the group test if the correlation coefficient distribution is different
% from 0 (see also the p_all output)
% D1 --> data matrix from list 1
% D2 --> data matrix for list 2 (if not empty)
% data_PET --> PET data matrix
% Resh --> Reshaped results matrix [file_index PET_map Correlation_result p-value file_path]


atlas_hdr = spm_vol(atlas);
atlas1 = spm_read_vols(atlas_hdr);
[a,b,c] = unique(atlas1(:));
a = a(a~=0);
atlas_vals = a(~isnan(a));

D1 = mean_time_course(list1,atlas, atlas_vals);

if ~isemptycell(list2)
    D2 = mean_time_course(list2,atlas, atlas_vals); 
else
    D2 = [];
end

data_PET = mean_time_course(files_PET,atlas,atlas_vals);

if options(4)==1 % adjust for structural correlation
    path_T1 = fullfile(fileparts(which('spm')),'tpm','TPM.nii,1');
    T1 =  mean_time_course({path_T1},atlas, atlas_vals); 
else
    T1 = '';
end

switch options(1)
    % opt_comp = 1 --> es between
    % opt_comp = 2 --> es within
    % opt_comp = 3 --> mean list 1
    % opt_comp = 4 --> list 1 each
    % opt_comp = 5 --> ind z-score list 1 to list 2
    % opt_comp = 6 --> pair-wise difference list 1 to list 2
    % opt_comp = 7 --> ind z-scores from list 1
    % opt_comp = 8 --> list 1 each compares against null distribution of
    % correlation coefficients
    case 1 % Cohen's d between groups
        m_D1 = mean(D1);
        std_D1 = std(D1);
        m_D2 = mean(D2);
        std_D2 = std(D2);
        data = (m_D1-m_D2)./sqrt((std_D1.^2+std_D2.^2)./2);
    case 2 % Cohen's d within group change
        delta_d = D1-D2;
        data = mean(delta_d)./std(delta_d);
    case 3 % mean list 1
        if size(D1,1)==1
            data = D1;
        else
            data = mean(D1);
        end
    case 4 % list 1 with PET data
        data = D1;
    case 5 % compute z-score list 1 relative to list 2
        m_D2 = mean(D2);
        std_D2 = std(D2);
        data = (D1 - repmat(m_D2,size(D1,1),1))./repmat(std_D2,size(D1,1),1);
    case 6 % pair-wise differences list 1 - list 2
        data = D1 - D2;
    case 7 % leave one out
        for i = 1:size(D1,1)
            data_i = D1(i,:);
            data_z = D1([1:i-1 i+1:end],:);
            data(i,:) = (data_i-mean(data_z))./std(data_z);
        end
    case 8
        data = D1;
end

% if options(2) == 3 % 1 = Spearman, 2 = Pearson, 3 = multiple linear regression
%     [res,p_all,stats] = correlateModalities(data,data_PET, options(2));
% else
if options(4)==1
    [res,p_all,stats] = correlateModalities(data,data_PET,options,T1);
else 
    [res,p_all,stats] = correlateModalities(data,data_PET,options);
end
% end
[Resh] = reshape_res(options,res,p_all, files_PET,list1);
    a = 1;
    if options(1)<4
        Y = zeros(size(atlas1));
        for i = 1:length(atlas_vals)
            Y(atlas1(:)==atlas_vals(i)) = data(i);
        end
        atlas_hdr.fname = image_save;
        atlas_hdr.pinfo(1) = 0.001;
        atlas_hdr.dt = [16 0];
        spm_write_vol(atlas_hdr,Y);
    end
end

function [res,p_all,stats] = correlateModalities(data,data_PET,opts,T1)
    

for i = 1:size(data,1)
    switch opts(2)
        case 1 % Spearman correlation
            for j = 1:size(data_PET,1)
                if opts(4)==1
                    [r,p] = partialcorr(data(i,:)',data_PET(j,:)',T1','type','Spearman');
                else
                    [r,p] = corr(data(i,:)',data_PET(j,:)','type','Spearman');
                end
                    corrOrig(i,j) = r;
                    p_ind(i,j) = p;
            end
            res_ind = fishers_r_to_z(corrOrig);
            stats.corrOrig = corrOrig;
            stats.p_ind = p_ind;
            stats.res_ind = res_ind;
            
            
        case 2 % Pearson correlation
            for j = 1:size(data_PET,1)
                if opts(4)==1
                    [r,p] = partialcorr(data(i,:)',data_PET(j,:)',T1','type','Pearson');
                else
                    [r,p] = corr(data(i,:)',data_PET(j,:)','type','Pearson');
                end
                    corrOrig(i,j) = r;
                    p_ind(i,j) = p;
            end 
            res_ind = fishers_r_to_z(corrOrig);
            stats.corrOrig = corrOrig;
            stats.p_ind = p_ind;
            stats.res_ind = res_ind;
        
        case 3 % multiple linear regresion
            y = data(i,:)';
            if opts(4)==1
                X = [data_PET' T1'];
                ind_PET = 2:size(X,2);
            else
                X = [data_PET'];
                ind_PET = 2:size(X,2)+1;
            end
            
            
            Stats = regstats(y,X);
            res_ind(i,:) = Stats.beta(ind_PET)';
            stats.res_ind(i,:) = Stats.beta(ind_PET)';
            stats.ind(i).p_ftest_full = Stats.fstat.pval;
            stats.ind(i).f_ftest_full = Stats.fstat.f;
            stats.ind(i).t_all = Stats.tstat.t;
            stats.p_ind(i,:) = Stats.tstat.pval(ind_PET)';

            


    end
end

    % opt_comp = 1 --> es between
    % opt_comp = 2 --> es within
    % opt_comp = 3 --> mean list 1
    % opt_comp = 4 --> list 1 each
    % opt_comp = 5 --> ind z-score list 1 to list 2
    % opt_comp = 6 --> pair-wise difference list 1 to list 2
    % opt_comp = 7 --> ind z-scores from list 1
    % opt_comp = 8 --> list 1 each compares against null distribution of
    % correlation coefficients
 switch opts(1)
     case {1,2,3,4} 
     res = res_ind;
     p_all  = stats.p_ind;
     
     case {5,6,7,8}
     res = mean(res_ind);
     [h,p_all,ci] = ttest(res_ind);
     stats.ci95 = ci;
 end

end

function [Resh]  = reshape_res(options,res,p_all, files_PET,list1)


for i = 1:length(files_PET);
   [path,name] = fileparts(files_PET{i});
   tt = regexp(name,'_','Split');
   Rec_list{i}=tt{1};
end

switch options(1)
    case 1
        ff = 'Effect Size';
    case 2
        ff = 'Effect Size (within subject)';
    case 3
        ff = 'Mean list';
    case 4
        ff = 'List 1 image';
    case 5
        ff = 'Z-score file';
    case 6
        ff = 'Delta';
    case 7
        ff = 'Z-score loo';
    case 8
        ff = 'List 1 all against null';
        
end

switch options(2)
    case 1
        if options(1)<4
            Resh = {ff 'PET Map' 'Fisher''s z (Spearman rho)' 'p-val (parametric)' 'File'};
        else
            Resh = {ff 'PET Map' 'Mean Fisher''s z (Spearman rho)' 'p-val (parametric)' 'File'};
        end
    case 2
        if options(1)<4
            Resh = {ff 'PET Map' 'Fisher''s z (Pearson r)' 'p-val (parametric)' 'File'};
        else
            Resh = {ff 'PET Map' 'Mean Fisher''s z (Pearson r)' 'p-val (parametric)' 'File'};
        end
    otherwise
        Resh = {'File' 'PET Map' 'Mean Beta' 'p-val (Regression analysis)' 'File'};
end

for i = 1:size(res,1)
    for j = 1:size(res,2)
       Resh{end+1,1} = ff; %[ff ' ' num2str(i)];
       Resh{end,2} = Rec_list{j};
       Resh{end,3} = res(i,j);
       if p_all(i,j)<.001
           p_i_j = '<.001';
       else
           p_i_j = num2str(p_all(i,j));
       end
       Resh{end,4} = p_i_j;
       Resh{end,5} = list1{i};
    end
end
end
