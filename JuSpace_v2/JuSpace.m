% App designed with MATLAB App Designer
classdef JuSpace < matlab.apps.AppBase

    % Components declaration
    properties (Access = public)
        UIFigure                       matlab.ui.Figure
        TabGroup                       matlab.ui.container.TabGroup
        Tab1                           matlab.ui.container.Tab
        Tab2                           matlab.ui.container.Tab
        Tab3                           matlab.ui.container.Tab
        LogoImage                      matlab.ui.control.Image

        % Tab 1 components
        AtlasLabel                     matlab.ui.control.Label
        AtlasDropDown                  matlab.ui.control.DropDown
        SelectAtlasButton              matlab.ui.control.Button
        AnalysisLabel                  matlab.ui.control.Label
        AnalysisDropDown               matlab.ui.control.DropDown
        FirstSetButton                 matlab.ui.control.Button
        SecondSetButton                matlab.ui.control.Button
        SelectSaveDirButton            matlab.ui.control.Button
        SaveDirLabel                   matlab.ui.control.Label
        FirstSetListBox                matlab.ui.control.ListBox
        SecondSetListBox               matlab.ui.control.ListBox
        FirstSetLabel                  matlab.ui.control.Label
        SecondSetLabel                 matlab.ui.control.Label
        AnalysisTooltipLabel           matlab.ui.control.Label
        NameSaveField                  matlab.ui.control.EditField
        NameSaveFieldLabel             matlab.ui.control.Label
        
        % Tab 2 components
        NeuroTemplateLabel             matlab.ui.control.Label
        NeuroTemplateListBox           matlab.ui.control.ListBox
        CellularMarkerLabel            matlab.ui.control.Label
        CellularMarkerListBox          matlab.ui.control.ListBox
        MetricsLabel                   matlab.ui.control.Label
        MetricsListBox                 matlab.ui.control.ListBox
        RunAnalysisButton              matlab.ui.control.Button
        PermutationsField              matlab.ui.control.EditField
        PermutationsFieldLabel         matlab.ui.control.Label
        T1CheckBox                     matlab.ui.control.CheckBox
        InstructionsLabel              matlab.ui.control.Label
      
        
        % Tab 3 components
        PlotSelected                   matlab.ui.control.Button
        FigureNeuro                    matlab.ui.control.UIAxes
        ResultsTable                   matlab.ui.control.Table
        SwitchScatterLabel1            matlab.ui.control.Label
        SwitchScatterLabel2            matlab.ui.control.Label
        SwitchScatter                  matlab.ui.control.Switch
        SaveFigureButton               matlab.ui.control.Button
        RankPlotCheckBox               matlab.ui.control.CheckBox
        ResultsLabel                   matlab.ui.control.Label
        xaxisListLabel                 matlab.ui.control.Label
        xaxisList                      matlab.ui.control.ListBox
        yaxisList                      matlab.ui.control.ListBox
        FDRCheckBox                    matlab.ui.control.CheckBox
        
        % Navigation buttons
        NextTabButton                  matlab.ui.control.Button
        PrevTabButton                  matlab.ui.control.Button
    end
    
    properties (Access = private)
        list_PET;
        list_cell;
        files_set1;
        files_set2;
        atlas;
        Results;
        Tab3Enabled = false;
        dir_tool;
    end
    
    
    methods (Access = private)

        function tabChanged(app, event)
            selectedTab = event.NewValue;

                if selectedTab == app.Tab3 && ~app.Tab3Enabled
                    % Block access: force back to Tab1
                    uialert(app.UIFigure, 'Please set-up and run the analysis first.', 'Tab Locked');
                    app.TabGroup.SelectedTab = app.Tab1;  % revert selection
                end
            check_inputs(app)
        end
        
        % Navigate tabs
        function nextTab(app, ~)  % <-- added 'event' here
            idx = find(app.TabGroup.Children == app.TabGroup.SelectedTab);
            if idx < 3
                app.TabGroup.SelectedTab = app.TabGroup.Children(idx + 1);
            end
            drawnow;
            check_inputs(app)
        end

        function prevTab(app, ~)
            idx = find(app.TabGroup.Children == app.TabGroup.SelectedTab);
            if idx > 1
                app.TabGroup.SelectedTab = app.TabGroup.Children(idx - 1);
            end
            check_inputs(app)
        end
        
        

        % Select custom Atlas
        function selectAtlas(app, ~, ~)
            [file, path] = uigetfile('*.nii;*.img', 'Select Atlas Image');
            atlas_all = app.atlas;
            atlas_all(end+1).name = file;
            atlas_all(end).folder = path;
            app.atlas = atlas_all;
          
            if isequal(file,0)
                app.AtlasDropDown.Value = 'm_labels_Neuromorphometrics.nii';
            else
                app.AtlasDropDown.Items = [app.AtlasDropDown.Items, {file}];
                app.AtlasDropDown.Value = file;
            end
            check_inputs(app)
        end
        
        %  Study design options
        function updateAnalysisVisibility(app)
        % update visibility of set 2
        
        study_design_opt = find(ismember(app.AnalysisDropDown.Items,app.AnalysisDropDown.Value))-1;
        set2_opt = [1,2,5,6];
              
                if ismember(study_design_opt, set2_opt)
                    app.SecondSetButton.Enable = 'on';
                else
                    app.SecondSetButton.Enable = 'off';
                    app.SecondSetListBox.Items = {''};
                    app.files_set2 = {''};
                    app.SecondSetButton.BackgroundColor = [0.5 0.7 0.9];
                end

            % Update tooltip based on dropdown selection
            
                switch study_design_opt
                    % opt_comp = 1 --> es between
                    % opt_comp = 2 --> es within
                    % opt_comp = 3 --> mean list 1
                    % opt_comp = 4 --> list 1 each
                    % opt_comp = 5 --> ind z-score list 1 to list 2
                    % opt_comp = 6 --> pair-wise difference list 1 to list 2
                    % opt_comp = 7 --> ind z-scores from list 1
                    % opt_comp = 8 --> list 1 each compares against null distribution of
                    % correlation coefficients
                    case 1
                        app.AnalysisTooltipLabel.Text = sprintf('Computes and uses Cohen''s d effect size per region for set 1 versus set 2');
                    case 2
                        app.AnalysisTooltipLabel.Text = sprintf('Computes pre-post Cohen''s d effect size \nper region for set 1 relative to set 2. \nNumber of selected images must be the same for set 1 and 2');
                    case 3
                        app.AnalysisTooltipLabel.Text = sprintf('Uses mean per region from set 1');
                    case 4
                        app.AnalysisTooltipLabel.Text = sprintf('Tests each image from set 1 against a null distribution');
                    case 5
                        app.AnalysisTooltipLabel.Text = sprintf('Converts each image in set 1 to z-scores relative to images in set 2');
                    case 6
                        app.AnalysisTooltipLabel.Text = sprintf('Computes pair-wise differences for set 1 versus set 2 \nand tests these difference maps against null distribution. \nNumber of selected images must be the same for set 1 and 2');
                    case 7
                        app.AnalysisTooltipLabel.Text = sprintf('Converts each image in set 1 to z-scores relative to \nall other images in set 1');
                    case 8  
                        app.AnalysisTooltipLabel.Text = sprintf('Tests the distribution of spatial correlation of all \nimages from set 1 against a null distribution');
                    otherwise
                        app.AnalysisTooltipLabel.Text = sprintf('Please select an analysis option');
                end
            check_inputs(app)
        end
        
        
                % Select First Image Set
                function selectFirstSet(app, ~, ~)
            
            app.files_set1 = cellstr(spm_select(Inf,'image','Select files for set 1'));
            for i = 1:length(app.files_set1)
                [~,file] = fileparts(app.files_set1{i});
                set1{i} = file;
            end
            
            app.UIFigure.Visible = 'off';
            app.UIFigure.Visible = 'on'; % keeps figure in the foreground
            if isequal(set1,0)
                app.FirstSetListBox.Items = {''};
            else
                if ischar(set1)
                    set1 = {set1};
                end
                app.FirstSetListBox.Items = set1;
                app.FirstSetListBox.Visible = 'on';
                app.FirstSetLabel.Visible = 'on';
                app.FirstSetButton.BackgroundColor = [0.4 0.8 0.6];
            end
            check_inputs(app)
        end

        % Select Second Image Set (optional)
        function selectSecondSet(app, ~, ~)
            
            ana_opt = find(ismember(app.AnalysisDropDown.Items,app.AnalysisDropDown.Value))-1;


            app.files_set2 = cellstr(spm_select(Inf,'image','Select files for set 1'));
            
            if ismember(ana_opt,[2 6]) && length(app.files_set1) ~= length(app.files_set2)
                uialert(app.UIFigure, 'For analysis options 2 and 6 the number of images in set 2 must match the number selected in set 1', 'Check selection');
            else
                for i = 1:length(app.files_set2)
                    [~,file] = fileparts(app.files_set2{i});
                    set2{i} = file;
                end
                app.UIFigure.Visible = 'off'; 
                app.UIFigure.Visible = 'on'; 
                if isequal(set2,0)
                    app.SecondSetListBox.Visible = 'off';
                    app.SecondSetLabel.Visible = 'off';
                else
                    if ischar(set2)
                        set2 = {set2};
                    end
                    app.SecondSetListBox.Items = set2;
                    app.SecondSetListBox.Visible = 'on';
                    app.SecondSetLabel.Visible = 'on';
                    app.SecondSetButton.BackgroundColor = [0.4 0.8 0.6];
                end
                check_inputs(app)
            end
        end
        
        % Select save directory
        function selectSaveDir(app, ~, ~)
            saveDir = uigetdir;
            if saveDir ~= 0
                app.SaveDirLabel.Text = saveDir;
                app.SelectSaveDirButton.BackgroundColor = [0.4 0.8 0.6];
            end
            app.UIFigure.Visible = 'off'; 
            app.UIFigure.Visible = 'on';   
            drawnow;
            check_inputs(app);
        
        end
        
        
        function run_analysis(app,~)
            
            d = uiprogressdlg(app.UIFigure, 'Title', 'Please Wait', 'Message', 'Running analysis...', 'Indeterminate', 'on', 'Cancelable', 'off');

            try
            % get selections from tab 1 (Inputs)

                try
                d.Message = 'Loading all settings ...';
                atlases = app.atlas;
                atlas_sel = app.AtlasDropDown.Value;
                atlas_all = app.AtlasDropDown.Items;
                [~,ind_atlas] = ismember(atlas_sel, atlas_all);
                atlas = fullfile(atlases(ind_atlas).folder,atlases(ind_atlas).name);
                ana_opt = find(ismember(app.AnalysisDropDown.Items,app.AnalysisDropDown.Value))-1;
                list1 = app.files_set1;
                list2 = app.files_set2;
                dir_save = app.SaveDirLabel.Text;
                name_save = app.NameSaveField.Value;
                app.SwitchScatter.Value = {'Bar plot'};
                app.xaxisList.Items = {''};
                app.yaxisList.Items = {''};
    
    
                % get selections from tab 2 (Templates and Metrics)
                [ind_PET,ind_cell,ind_ana] = SelectTemplates(app);
                files_PET = app.list_PET(ind_PET);
                files_cell = app.list_cell(ind_cell);
                files_PET = [files_PET;files_cell];
                opt_perm  = 0;
                opt_spat_perm = 0;
                Nperm = str2num(app.PermutationsField.Value);
                
    
    
                switch ana_opt
                    case {1,2,5,6}
                        opt_perm = 1;
                    case {3,4,7,8}
                        opt_spat_perm = 1;
                end
                image_save = fullfile(dir_save, [name_save '.nii']);
                
                opt_T1 = app.T1CheckBox.Value;
               
                options = [ana_opt ind_ana opt_perm opt_T1 opt_spat_perm];
                
                catch ME
                    uialert(app.UIFigure, getReport(ME), 'Analysis Error');
                    fid = fopen('error_log.txt', 'a');
                    fprintf(fid, '[%s] ERROR: %s\n', char(datetime), ME.message);
                    fprintf(fid, '%s\n\n', getReport(ME, 'extended'));
                    fclose(fid);
                end
                d.Message = 'Computing spatial correlations...';
                try
                    Results = compute_DomainGauges(list1,list2,files_PET,atlas, options,image_save);
                catch ME
                    fid = fopen('error_log.txt', 'a');
                    fprintf(fid, '[%s] ERROR: %s\n', char(datetime), ME.message);
                    fprintf(fid, '%s\n\n', getReport(ME, 'extended'));
                    fclose(fid);
                end
    
                opt_for_perm = [1,2,5,6];
                opt_for_spat_perm = [3, 4, 7, 8];

                
                try
                    if options(3)==1 && ismember(options(1),opt_for_perm)% && options(2)~=3
                        d.Message = 'Computing exact p-value...';
                        disp('Computing exact p-value');
                       [p_exact,~] = compute_exact_pvalue(Results.data_set1,Results.data_set2,Results.data_PET,Results.res,Nperm,options,Results.T1,Results.stats);
                       Results.Resh(:,end+1) = [{'p_exact'}; num2cell_my(p_exact')];
                    end
        
                    if options(5)==1 && ismember(options(1),opt_for_spat_perm)
                        d.Message = 'Computing exact spatial p-value. This option may take quite a while if the null maps are not yet precomputed';
                        disp('Computing exact spatial p-value')
                        [p_exact,~,~] = compute_exact_spatial_pvalue(Results.data_set1,Results.data_PET,atlas,Results.res,Nperm,options,files_PET, Results.T1,Results.stats,d);
                        Results.Resh(:,end+1) = [{'p_exact_spatial'}; num2cell_my(p_exact')];
                    end
                    [h_sig, crit_p, adj_p_fdrBH] = fdr_bh(p_exact, 0.05);
                    Results.Resh(:,end+1) = [{'p_fdr_BH'}; num2cell_my(adj_p_fdrBH)]; 
                catch ME
                    fid = fopen('error_log.txt', 'a');
                    fprintf(fid, '[%s] ERROR: %s\n', char(datetime), ME.message);
                    fprintf(fid, '%s\n\n', getReport(ME, 'extended'));
                    fclose(fid);
                end
                d.Message = 'Saving and visualizing results...';
                if options(2)<3
                    Results.Resh = Results.Resh(:,[1 6 7 2 3 4 5]);
                else
                    Results.Resh = Results.Resh(:,[1 5 6 2 3 4]);
                end
                Results.p_exact = p_exact;
                Results.p_exact_fdr_BH = adj_p_fdrBH;
                Results.dir_save = dir_save;
                Results.Nperm = Nperm;
                Results.options = options;
                Results.atlas = atlas;
                Results.filesPET = files_PET;
                Results.set1_images = list1;
                Results.set2_images = list2;
                app.Results = Results;
                app.ResultsTable.ColumnName = Results.Resh(1,:);
                
                data = Results.Resh(2:end,:); 
                app.ResultsTable.Data = data;
                save_results(Results,app,name_save);
                try
                    bar_plot(app,h_sig);
                    app.TabGroup.SelectedTab = app.TabGroup.Children(3);
                    app.Tab3Enabled = true;
                catch ME
                    fid = fopen('error_log.txt', 'a');
                    fprintf(fid, '[%s] ERROR: %s\n', char(datetime), ME.message);
                    fprintf(fid, '%s\n\n', getReport(ME, 'extended'));
                    fclose(fid);
                end

            catch ME2
                uialert(app.UIFigure, getReport(ME2), 'Analysis Error');
            end
            close(d);
        end
        
       function listBoxItemClicked(app, event)
            check_inputs(app);
       end
       
       function switchscatterbox(app)
           switch_plot =  app.SwitchScatter.Value;
           if strcmp(switch_plot,'Scatter plot')
               cla(app.FigureNeuro, 'reset'); 
               app.xaxisList.Enable = 'on';
               app.yaxisList.Enable = 'on';
               app.RankPlotCheckBox.Enable = 'on';
               app.FDRCheckBox.Enable = 'off';
               app.PlotSelected.Enable = 'on';
               all = app.Results.filesPET;
               for i = 1:length(all)
                   [~,name] = fileparts(all{i});
                   xlist_all{i,1} = name;
               end
               app.xaxisList.Items = xlist_all;
               if size(app.Results.data,1) > 1
                   app.yaxisList.Items = app.files_set1;
               else
                   app.yaxisList.Items = {'data modality'};
               end
           else
               app.FDRCheckBox.Enable = 'on';
               app.FDRCheckBox.Value = false;
               cla(app.FigureNeuro, 'reset'); 
               app.xaxisList.Enable = 'off';
               app.yaxisList.Enable = 'off';
               app.RankPlotCheckBox.Enable = 'off';
               app.PlotSelected.Enable = 'off';
               bar_plot(app);

           end
       end

       function addSignficanttoBoxPlot(app)

            plot_sig = app.FDRCheckBox.Value;
         
            if plot_sig
                try
                    h_sig = app.Results.p_exact_fdr_BH<0.05;
                    bar_plot(app,h_sig);
                catch
                    uialert(app.UIFigure, 'FDR significance not computed','Information');
                end 
            else
                bar_plot(app);
            end
       end
       
       function plotSelected(app)
           scatter_plot(app);
       end
    
        function exportWithPrint(app)
            % Create a new figure and axes
            
            plot_opt = app.SwitchScatter.Value;
            
            if strcmp(plot_opt,'Scatter plot')
                f = figure('Position',[400 400 800 600],'Visible', 'off');
            else
                if length(app.Results.filesPET)<3
                    f = figure('Position',[400 400 500 500],'Visible', 'off');
                elseif length(app.Results.filesPET)<10
                    f = figure('Position',[400 200 length(app.Results.filesPET).*150 700],'Visible', 'off');
                else
                    f = figure('Position',[400 200 length(app.Results.filesPET).*100 700],'Visible', 'off');
                end
            end

            ax = axes(f);

            % Copy content from App Designer UIAxes
            copyobj(allchild(app.FigureNeuro), ax);

            % Copy axis labels and limits
            ax.XLim = app.FigureNeuro.XLim;
            ax.YLim = app.FigureNeuro.YLim;
            xlabel(ax, app.FigureNeuro.XLabel.String);
            ylabel(ax, app.FigureNeuro.YLabel.String);
            title(ax, app.FigureNeuro.Title.String);

            set(ax, 'XTick', app.FigureNeuro.XTick);
            set(ax, 'XTickLabel', app.FigureNeuro.XTickLabel);
            set(ax, 'FontSize', app.FigureNeuro.FontSize);
            set(ax, 'XTickLabelRotation', app.FigureNeuro.XTickLabelRotation);
            try
            legend(ax, app.FigureNeuro.Legend.String, 'Location', 'eastoutside');
            catch
            end
            grid on
            set(gcf,'color','w');
            time_now = datestr(datetime('now'),'ddmmmyyyy_HHMMSS');
            
            if strcmp(plot_opt,'Scatter plot')
                print(f,fullfile(app.Results.dir_save,['Scatter_' app.NameSaveField.Value  '_' time_now '.png']),'-dpng','-r300');
            else
                print(f,fullfile(app.Results.dir_save,['Bar_' app.NameSaveField.Value  '_' time_now '.png']),'-dpng','-r300');
            end

            close(f);
        end
        
        function loadAnalysis(app)

            d = uiprogressdlg(app.UIFigure, 'Title', 'Please Wait', 'Message', 'Loading analysis...', 'Indeterminate', 'on', 'Cancelable', 'off');

            [file, path] = uigetfile('*.mat', 'Load Analysis');
            
            if isequal(file, 0)
                return; 
            end
            loaded = load(fullfile(path, file));
            
            if ~isfield(loaded, 'app_save')
                uialert(app.UIFigure, 'Invalid file format.', 'Error');
            return;
            end
            analysisData = loaded.app_save;
            
            % Restore Tab 1
            app.AtlasDropDown.Items = analysisData.AtlasDropDown.Items;
            app.AtlasDropDown.Value = analysisData.AtlasDropDown.Value;
            
            app.AnalysisDropDown.Items = analysisData.AnalysisDropDown.Items;
            app.AnalysisDropDown.Value = analysisData.AnalysisDropDown.Value;
            
            app.FirstSetListBox.Items = analysisData.FirstSetListBox.Items;
            app.FirstSetListBox.Value = analysisData.FirstSetListBox.Value;            
            app.SecondSetListBox.Items = analysisData.SecondSetListBox.Items;
            app.SecondSetListBox.Value = analysisData.SecondSetListBox.Value;  
            app.AnalysisTooltipLabel.Text = analysisData.AnalysisTooltipLabel.Text;
            app.NameSaveField.Value = analysisData.NameSaveField.Value;
            app.SaveDirLabel.Text = analysisData.SaveDirLabel.Text;
            app.files_set1 = analysisData.files_set1;
            app.files_set2 = analysisData.files_set2;
            app.list_PET = analysisData.list_PET;
            app.list_cell = analysisData.list_cell;
            app.atlas = analysisData.atlas;
            
            % Restore Tab 2
            app.NeuroTemplateListBox.Items = analysisData.NeuroTemplateListBox.Items;
            app.NeuroTemplateListBox.Value = analysisData.NeuroTemplateListBox.Value;
            app.CellularMarkerListBox.Items = analysisData.CellularMarkerListBox.Items;
            app.CellularMarkerListBox.Value = analysisData.CellularMarkerListBox.Value;
            app.MetricsListBox.Items = analysisData.MetricsListBox.Items;
            app.MetricsListBox.Value = analysisData.MetricsListBox.Value;   
            app.PermutationsField.Value = analysisData.PermutationsField.Value;   
            app.T1CheckBox.Value = analysisData.T1CheckBox.Value;               

            % Restore Tab 3
            if ~isempty(analysisData.ResultsTable.Data)
                cla(app.FigureNeuro, 'reset'); 
                app.TabGroup.SelectedTab = app.TabGroup.Children(3);
                app.Tab3Enabled = true;
                app.ResultsTable.ColumnName = analysisData.ResultsTable.ColumnName;
                app.ResultsTable.Data = analysisData.ResultsTable.Data;
                app.SwitchScatter.Value = analysisData.SwitchScatter.Value;  
                app.RankPlotCheckBox.Value = analysisData.RankPlotCheckBox.Value;  
                app.xaxisList.Items = analysisData.xaxisList.Items;
                app.xaxisList.Value = analysisData.xaxisList.Value;            
                app.yaxisList.Items = analysisData.yaxisList.Items;
                app.yaxisList.Value = analysisData.yaxisList.Value; 
                app.Results = loaded.Results;
                app.xaxisList.Enable = 'off';
                app.yaxisList.Enable = 'off';
                app.RankPlotCheckBox.Enable = 'off';
                app.PlotSelected.Enable = 'off';
                app.TabGroup.SelectedTab = app.TabGroup.Children(3);
                bar_plot(app);
            else
                app.TabGroup.SelectedTab = app.TabGroup.Children(1);
                app.Tab3Enabled = false;
            end
            close(d);
            drawnow;
            check_inputs(app);
            uialert(app.UIFigure, 'Analysis loaded.', 'Loaded');

        end
        
        function saveAnalysis(app)
         
            Results = app.Results;
            % Results.JuSpace_version = 'v2.0';
            % file_save = fullfile(path,file);
            % save(file_save,'Results','app_save');
            save_results(Results,app);
            uialert(app.UIFigure, 'Analysis saved.', 'Success');
        end
        

    end

    methods (Access = private)
        function createComponents(app)
            warning off
            
            if isdeployed
                 [~, ~] = system('path');
                 app.dir_tool = pwd; 
            else
                 app.dir_tool= fileparts(which('JuSpace'));
            end


     
            % Main UI figure
            app.UIFigure = uifigure('Name', 'JuSpace 2.0', 'Position', [100, 100, 900, 600], 'Color', [0.96 0.96 0.98], 'AutoResizeChildren','on');
            fileMenu = uimenu(app.UIFigure, 'Text', 'File');
            uimenu(fileMenu, 'Text', 'Load Existing Analysis...','MenuSelectedFcn', @(src, event) loadAnalysis(app));
            uimenu(fileMenu, 'Text', 'Save Analysis Settings...', 'MenuSelectedFcn', @(src, event) saveAnalysis(app));
            app.LogoImage = uiimage(app.UIFigure, 'ImageSource', fullfile(app.dir_tool,'splash.png'),'Position', [420, 1, 60, 60]);  

            % Tab Group
            app.TabGroup = uitabgroup(app.UIFigure, 'Position', [10 60 880 520],'SelectionChangedFcn', @(src,event) tabChanged(app, event));

            % Tab 1
            app.Tab1 = uitab(app.TabGroup, 'Title', 'Inputs');

            
            % Update AtlasDropDown items, keeping existing entries
        

            % Atlas options
            app.AtlasLabel = uilabel(app.Tab1, 'Text', 'Select Atlas:', 'Position', [20, 440, 100, 22], 'FontWeight', 'bold');
            app.AtlasDropDown = uidropdown(app.Tab1, 'Items', {'Atlas1'}, 'Position', [130, 440, 250, 22]);
            app.SelectAtlasButton = uibutton(app.Tab1, 'push', 'Text', 'Custom Atlas', 'Position', [400, 440, 120, 22], 'BackgroundColor',[0.5 0.7 0.9],'FontColor','white','ButtonPushedFcn', @app.selectAtlas);
            
            
            
            %Study design options
            study_designs = {'Select an option', '1) Effect size between groups', '2) Effect size within group','3) Mean from set 1', '4) Set 1 each image', '5) Individual z-scores for set 1 relative to set 2', '6) pair-wise difference set 1 relative to set 2','7) Leave-one-out from set 1', '8) Set 1 each compares against null distribution'}; 
            app.AnalysisTooltipLabel = uilabel(app.Tab1, 'Text', 'Description of the selected spatial correlation approach', 'Position', [400, 390, 400, 44],'FontAngle', 'italic', 'FontColor', [0.3 0.3 0.3]);
            app.AnalysisLabel = uilabel(app.Tab1, 'Text', 'Study Design:', 'Position', [20, 400, 100, 22], 'FontWeight', 'bold');
            app.AnalysisDropDown = uidropdown(app.Tab1, 'Items', study_designs, 'Position', [130, 400, 250, 22],'ValueChangedFcn', @(src, event) updateAnalysisVisibility(app));

            

            

            app.SelectSaveDirButton = uibutton(app.Tab1, 'push', 'Text', 'Select Save Directory', 'Position', [20, 350 360, 22],'BackgroundColor',[0.8 0.8 0.5],'FontColor','white' ,'ButtonPushedFcn',  @app.selectSaveDir);
            app.SaveDirLabel = uilabel(app.Tab1, 'Text', '', 'Position', [20, 320, 600, 22]);
            % Label for editable field
            app.NameSaveFieldLabel = uilabel(app.Tab1, 'Text', 'Name save', 'Position', [400, 370, 100, 22], 'FontWeight', 'bold');
            app.NameSaveField = uieditfield(app.Tab1, 'text', 'Position', [400, 350 200, 22], 'Value', ''); % initial value
            
            app.FirstSetButton = uibutton(app.Tab1, 'push', 'Text', 'Select First Set', 'Position', [20, 290, 120, 22],'BackgroundColor',[0.5 0.7 0.9],'FontColor','white','ButtonPushedFcn', @app.selectFirstSet);
            app.FirstSetLabel = uilabel(app.Tab1, 'Text', 'First Set Images:', 'Position', [20, 260, 340, 20],'FontWeight', 'bold');
            app.FirstSetListBox = uilistbox(app.Tab1,'Items',{''}, 'Position', [20, 30, 380, 220]);
            
            
            app.SecondSetButton = uibutton(app.Tab1, 'push', 'Text', 'Select Second Set', 'Position', [450, 290, 120, 22],'BackgroundColor',[0.5 0.7 0.9],'FontColor','white', 'ButtonPushedFcn', @app.selectSecondSet);
            app.SecondSetLabel = uilabel(app.Tab1, 'Text', 'Second Set Images:', 'Position', [450, 260, 340, 20],'FontWeight', 'bold');
            app.SecondSetListBox = uilistbox(app.Tab1,'Items',{''}, 'Position', [450, 30, 380, 220]);
            app.SecondSetButton.Enable = 'off';

            % Components in Tab 2 
            % Tab 2
            app.Tab2 = uitab(app.TabGroup, 'Title', 'Templates & Metrics');


            dir_PET = fullfile(app.dir_tool,'PETatlas');
            files_PET = select_con_maps_forfMRI_my(dir_PET,'PETatlas','.*.nii');
            
                for i = 1:length(files_PET)
                    [~,file] = fileparts(files_PET{i});
                    list_PET{i} = file;
                end
                
            app.list_PET = files_PET;
            app.NeuroTemplateLabel = uilabel(app.Tab2,'Text','PET/Neurotransmitter Templates:','Position',[20,470,200,22],'FontWeight','bold');
            app.NeuroTemplateListBox = uilistbox(app.Tab2,'Items',list_PET,'Position',[20,110,260,360],'Multiselect','on','ValueChangedFcn', @(src,event) listBoxItemClicked(app, event)); 
            app.NeuroTemplateListBox.Value = cell(0,0);
            
            dir_cell = fullfile(app.dir_tool,'CELLatlas');
            files_cell = select_con_maps_forfMRI_my(dir_cell,'CELLatlas','.*.nii');
                for i = 1:length(files_cell)
                    [path,file] = fileparts(files_cell{i});
                    list_cell{i} = file;
                end
            app.list_cell = files_cell;
            app.CellularMarkerLabel = uilabel(app.Tab2,'Text','Cellular Markers:','Position',[310,470,300,22],'FontWeight','bold');
            app.CellularMarkerListBox = uilistbox(app.Tab2,'Items',list_cell,'Position',[310,110,260,360],'Multiselect','on','ValueChangedFcn', @(src,event) listBoxItemClicked(app, event)); 
            app.CellularMarkerListBox.Value = cell(0,0);

            app.MetricsLabel = uilabel(app.Tab2,'Text','Select Analysis Type:','Position',[600,470,200,22],'FontWeight','bold');
            app.MetricsListBox = uilistbox(app.Tab2,'Items',{'Spearman','Pearson','Multiple Linear Regression'},'Position',[600,110,260,360],'Multiselect','off','ValueChangedFcn', @(src,event) listBoxItemClicked(app, event));
            app.MetricsListBox.Value = cell(0,0);
            
            
            app.PermutationsFieldLabel = uilabel(app.Tab2, 'Text', 'Number of permutations', 'Position', [20, 90, 200, 22], 'FontWeight', 'bold');
            app.PermutationsField = uieditfield(app.Tab2, 'text', 'Position', [20, 70, 150, 22], 'Value', '1000'); % initial value
            
            app.T1CheckBox = uicheckbox(app.Tab2, 'Text', 'Partial volume effect correction using T1 probabilities', 'Position', [200, 60, 350, 30]);
            
            app.RunAnalysisButton = uibutton(app.Tab2,'push','Text','Run Analysis','Position',[20,20,840,40],'BackgroundColor',[0.2 0.6 0.8],'FontColor','white', 'ButtonPushedFcn', @(src, event) run_analysis(app));
            app.RunAnalysisButton.Enable = 'off';
            
            app.InstructionsLabel = uilabel(app.Tab2, 'Text', sprintf('To (un)select multiple images \nkeep CTRL key pressed'), 'Position', [600, 70, 200, 33]);
            
            % Tab 3
            app.Tab3 = uitab(app.TabGroup, 'Title', 'Results');
            
            app.FigureNeuro = uiaxes(app.Tab3,'Position', [50, 60, 430, 430]);
            title(app.FigureNeuro, 'Results');
            app.ResultsLabel = uilabel(app.Tab3, 'Text', 'Results table', 'Position', [650, 470, 100, 22],'FontWeight','bold');
            app.ResultsTable = uitable(app.Tab3,'Position', [650, 20, 200, 450], 'Data', {},'ColumnName', {'Column 1', 'Column 2'},'RowName', {});
            
            % plot options
            app.SwitchScatter = uiswitch(app.Tab3, 'slider','Position', [60, 10, 45, 20],'Items', {'Bar plot', 'Scatter plot'}, 'ValueChangedFcn', @(src,event) switchscatterbox(app));
            app.RankPlotCheckBox = uicheckbox(app.Tab3, 'Text', 'Rank Plot', 'Position', [200, 8, 90, 20]);
            app.FDRCheckBox = uicheckbox(app.Tab3, 'Text', 'FDR sign.', 'Position', [290, 8, 100, 20],'ValueChangedFcn', @(src,event) addSignficanttoBoxPlot(app));

            app.SaveFigureButton  = uibutton(app.Tab3,'push','Text','Save figure','Position',[400,10,90,20], 'BackgroundColor',[0.2 0.6 0.8],'FontColor','white', 'ButtonPushedFcn',  @(src, event) exportWithPrint(app));
            
            app.xaxisListLabel = uilabel(app.Tab3, 'Text', 'Select plot items', 'Position', [530, 470, 100, 22],'FontWeight','bold');
            app.xaxisList = uilistbox(app.Tab3,'Items',{'x axis'},'Position',[530,260,110,210],'Multiselect','on'); 
            app.yaxisList = uilistbox(app.Tab3,'Items',{'y axis'},'Position',[530,40,110,210],'Multiselect','on'); 
            
            app.PlotSelected  = uibutton(app.Tab3,'push','Text','Plot selected','Position',[530,20,110,20], 'ButtonPushedFcn',  @(src, event) plotSelected(app));
            app.xaxisList.Enable = 'off';
            app.RankPlotCheckBox.Enable = 'off';
            app.PlotSelected.Enable = 'off';
            app.yaxisList.Enable = 'off';
            
            % Navigation Buttons
            app.PrevTabButton = uibutton(app.UIFigure, 'push', 'Text', 'Previous', 'Position', [10, 20, 100, 30], 'ButtonPushedFcn', @(src,event) app.prevTab(event));
            app.NextTabButton = uibutton(app.UIFigure, 'push', 'Text', 'Next', 'Position', [790, 20, 100, 30], 'ButtonPushedFcn', @(src,event) app.nextTab(event));
            
        end
    end
        

    methods (Access = public)
        function app = JuSpace
            createComponents(app);
            load_atlases(app);
        end
    end

end

    function load_atlases(app)  % checks for available atlases and sets the default
        dir_juspace = app.dir_tool;
        % dir_juspace = fileparts(which('JuSpace')); 
        dir_atlas = dir(fullfile(dir_juspace, 'atlas','*.nii')); % Adjust file extension if needed
        app.atlas = dir_atlas;
        atlases = {dir_atlas.name};
        app.AtlasDropDown.Items = atlases;
        app.AtlasDropDown.Value = 'm_labels_Neuromorphometrics.nii';
    end
    
    
    function [check_analysis] = check_inputs(app)
            %check inputs
            % tab 1
              study_design = find(ismember(app.AnalysisDropDown.Items,app.AnalysisDropDown.Value))-1;
              study_design_opt = study_design>0;
              set2_opt = [1,2,5,6];
              list1_opt = sum(~isemptycell(app.FirstSetListBox.Items))>0;

              save_dir_opt = isdir(app.SaveDirLabel.Text);
              check_all = [study_design_opt list1_opt save_dir_opt];

              if ismember(study_design, set2_opt)
                list2_opt = sum(~isemptycell(app.SecondSetListBox.Items))>0;
                check_all = [check_all list2_opt];
              end

              [ind1,ind2,ind_ana] = SelectTemplates(app);
              
              ind_check = ~isempty([ind1 ind2]);
              ind_ana_check = ~isempty(ind_ana);
              
              check_all = [check_all ind_check ind_ana_check];
              
              if true(check_all)
                  app.RunAnalysisButton.Enable = 'on';
                  drawnow;
              else
                  app.RunAnalysisButton.Enable = 'off';
                  drawnow;
              end

    end
    
    function  [ind_PET,ind_cell,ind_ana] = SelectTemplates(app)
                selectedItems = app.NeuroTemplateListBox.Value;
                allItems = app.NeuroTemplateListBox.Items;
                [~, ind_PET] = ismember(selectedItems, allItems);

                selectedItemsCell = app.CellularMarkerListBox.Value;
                allItemsCell = app.CellularMarkerListBox.Items;
                [~, ind_cell] = ismember(selectedItemsCell, allItemsCell);
                
                selectedItems = app.MetricsListBox.Value;
                allItems = app.MetricsListBox.Items;
                [~, ind_ana] = ismember(selectedItems, allItems);
    end
    
    
    function save_results(Results,app,name_save)
       
            Results.JuSpace_version = 'v2.0';

            % Save Tab 1
            app_save.AtlasDropDown.Items = app.AtlasDropDown.Items;
            app_save.AtlasDropDown.Value = app.AtlasDropDown.Value;
            app_save.AnalysisDropDown.Items = app.AnalysisDropDown.Items;
            app_save.AnalysisDropDown.Value = app.AnalysisDropDown.Value;
            app_save.FirstSetListBox.Items = app.FirstSetListBox.Items;
            app_save.FirstSetListBox.Value = app.FirstSetListBox.Value;            
            app_save.SecondSetListBox.Items = app.SecondSetListBox.Items;
            app_save.SecondSetListBox.Value = app.SecondSetListBox.Value;  
            app_save.AnalysisTooltipLabel.Text = app.AnalysisTooltipLabel.Text;
            app_save.NameSaveField.Value = app.NameSaveField.Value;
            app_save.SaveDirLabel.Text = app.SaveDirLabel.Text;
            app_save.files_set1 = app.files_set1;
            app_save.files_set2 = app.files_set2;
            app_save.list_PET = app.list_PET;
            app_save.list_cell = app.list_cell;
            app_save.atlas = app.atlas;
            
                        
            % Save Tab 2
            try
            app_save.NeuroTemplateListBox.Items = app.NeuroTemplateListBox.Items;
            app_save.NeuroTemplateListBox.Value = app.NeuroTemplateListBox.Value;
            app_save.CellularMarkerListBox.Items = app.CellularMarkerListBox.Items;
            app_save.CellularMarkerListBox.Value = app.CellularMarkerListBox.Value;
            app_save.MetricsListBox.Items = app.MetricsListBox.Items;
            app_save.MetricsListBox.Value = app.MetricsListBox.Value;   
            app_save.PermutationsField.Value = app.PermutationsField.Value;   
            app_save.T1CheckBox.Value = app.T1CheckBox.Value;  
            catch
            end
            time_now = datestr(datetime('now'),'ddmmmyyyy_HHMMSS');
            % Save Tab 3
            try
            app_save.TabGroup.SelectedTab = app.TabGroup.Children(3);
            app_save.Tab3Enabled = true;
            app_save.ResultsTable.ColumnName = app.Results.Resh(1,:);
            app_save.ResultsTable.Data = app.ResultsTable.Data;
            app_save.SwitchScatter.Value = app.SwitchScatter.Value;  
            app_save.RankPlotCheckBox.Value = app.RankPlotCheckBox.Value;  
            app_save.xaxisList.Items = app.xaxisList.Items;
            app_save.xaxisList.Value = app.xaxisList.Value;            
            app_save.yaxisList.Items = app.yaxisList.Items;
            app_save.yaxisList.Value = app.yaxisList.Value; 
            app_save.Results = app.Results;
            app_save.xaxisList.Enable = 'off';
            app_save.yaxisList.Enable = 'off';
            app_save.RankPlotCheckBox.Enable = 'off';
            app_save.PlotSelected.Enable = 'off';
           
            writecell(app.Results.Resh,fullfile(Results.dir_save,['ResultsTable_' name_save '_' time_now '.csv']),'Delimiter',';');
            catch
            end

            if ~exist('name_save','var') 
                [file, path] = uiputfile('*.mat', 'Save Analysis to...');
                    if isequal(file, 0)
                        return;  % User canceled
                    end
                save(fullfile(path,file),'Results','app_save');
            else
                name_save = add_suffix(name_save,Results.options);
                file_save = fullfile(Results.dir_save,[name_save '.mat']);
                save(file_save,'Results','app_save');
            end
    end
    

    
    function bar_plot(app,h_sig)
        cla(app.FigureNeuro, 'reset'); 
        Results = app.Results;
        name_save = app.NameSaveField.Value;
        res = Results.stats.res_ind;
        ind_plot =[];
        for i = 1:size(res,1)
            for j = 1:size(res,2)
               ind_plot(end+1,1) = j;
            end
        end

        vals_plot= res; %cell2num_my(Resh(2:end,3));
        for i = 1:length(Results.filesPET)
           all_i = vals_plot(:,i);%vals_plot(ind_plot==i);

           size_y_xi = sum(ind_plot==i);
           y_dist = abs(vals_plot(ind_plot==i) - mean(vals_plot(ind_plot==i)));
           y_dist_inv = abs(y_dist -max(y_dist))./max(abs(y_dist -max(y_dist)));
           if isnan(y_dist_inv)
               y_dist_inv = rand(size(y_dist_inv));
           end

           x_n_n(ind_plot==i) = ind_plot(ind_plot==i) + 0.7.*(rand(size_y_xi,1)-0.5).*y_dist_inv;
           all_ii = all_i(~isinf(all_i));
           m_x(i) = mean(all_ii);
           std_x(i,1) = 1.95996.*std(all_ii)./sqrt(length(all_ii));
        end
        
        
        % plot specs

        n_files = length(Results.set1_images);
        marker_color = [0 0.7 1];
        opt_plot = 1;
        if opt_plot == 1
            bar_color = generate_colors_nice_my(length(Results.filesPET));
        else
            bar_color = generate_colors_blue_my(length(Results.filesPET));
        end

        opt_comp = Results.options(1);
        if opt_comp<3 || size(Results.data,1)==1
            h2 = bar(app.FigureNeuro,1:length(Results.filesPET),diag(m_x),0.9,'stacked','EdgeColor','none');
            if min(m_x)>0
                min_m_x = 0;
            else
                min_m_x = min(m_x);
            end

            if max(m_x)<0
                max_m_x = 0;
            else
                max_m_x = max(m_x);
            end
            ylim(app.FigureNeuro,[min_m_x.*1.1 max_m_x.*1.1]);

        else

            h2 = boxplot(app.FigureNeuro, vals_plot,'widths',0.9,'outliersize',2);
            set(app.FigureNeuro,'Visible','on');
        %     set(gca,'XTickLabel',Rec_list);
            h3 = findobj(app.FigureNeuro,'Tag','Box');
            flip_color = flip(bar_color);
            for j=1:length(h3)
                patch(app.FigureNeuro,get(h3(j),'XData'),get(h3(j),'YData'),flip_color(j,:),'EdgeColor','none');
            end
            hold(app.FigureNeuro, 'on');
            h2 = boxplot(app.FigureNeuro, vals_plot,'widths',0.9,'outliersize',2,'Colors',zeros(3,3));
            set(h2(7,:),'Visible','off');
            lines = findobj(app.FigureNeuro, 'type', 'line', 'Tag', 'Median');
            set(lines, 'Color', 'black');
            grid(app.FigureNeuro,'on');
        end

        Rec_list = Results.Resh(2:end,1);
        app.FigureNeuro.XTick = 1:length(Rec_list);
        
        if length(Results.filesPET)>1
            app.FigureNeuro.XTickLabel = Rec_list;
            app.FigureNeuro.FontSize = 16;
            app.FigureNeuro.XTickLabelRotation = 90;
        else
            app.FigureNeuro.XTickLabel = Rec_list;
            app.FigureNeuro.FontSize = 16;
        end

         vv = vals_plot';
         vv2 = vv(:);
        hold(app.FigureNeuro, 'on');
        for i = 1:length(Results.filesPET)
            if size(Results.data,1)==1
                set(h2(i),'facecolor',bar_color(i,:),'edgecolor','none');
            end
            if opt_comp>=4

                if n_files>1
                plot(app.FigureNeuro, x_n_n(ind_plot==i),vv2(ind_plot==i),'o','MarkerSize',8,'MarkerFaceColor',bar_color(i,:),'MarkerEdgeColor','black','LineWidth',1);
                end
            end
        end

        if opt_comp~=4 && opt_comp ~= 7 
            if exist('h_sig','var') 
                ind_sig = find(h_sig);
                ind_sig = ind_sig(ind_sig<=max(app.FigureNeuro.XTick));

                if ~exist('max_m_x','var')
                    max_m_x = max(vals_plot(:));
                    min_m_x = min(vals_plot(:));
                end
                % Add asterisks+
                    for i = 1:length(ind_sig)
                            text(app.FigureNeuro,ind_sig(i), max_m_x+0.04, '*', ...
                                'HorizontalAlignment', 'center',...
                                'FontSize', 40, 'FontWeight', 'bold');
                    end
                ylim(app.FigureNeuro,[min_m_x.*1.1 max_m_x+0.20]);
            end 
        end

       app.FigureNeuro.YLimMode = 'manual';

       opt_ana = Results.options(2);
        switch opt_ana
            case 1
                ylabel(app.FigureNeuro,'Fisher''s z (Spearman rho)','fontweight','bold');
            case 2
                ylabel(app.FigureNeuro,'Fisher''s z (Pearson r)','fontweight','bold');
            case 3
                ylabel(app.FigureNeuro,'Standardized beta coeffiecient','fontweight','bold');
        end
        app.FigureNeuro.Color = 'white';
        hold(app.FigureNeuro, 'off');

    end
    
    function scatter_plot(app)

        cla(app.FigureNeuro, 'reset'); 
        opt_ana = app.Results.options(1);
        data_list = app.Results.set1_images;
        data_all = app.Results.data;
        All_PET = app.xaxisList.Items;
        
        Selected_PET = app.xaxisList.Value;
        [~, ind_sel] = ismember(Selected_PET, All_PET);
        
        if size(data_all,1)>1
             selected_data = app.yaxisList.Value;
            ind_sel_data = find(ismember(data_list, selected_data));
        else
            ind_sel_data = 1;
        end
        
        Rec_list_all = app.Results.Resh(2:end,:);
        Rec_list = Rec_list_all(ind_sel);
        data_PET = app.Results.data_PET(ind_sel,:);
        colors = generate_colors_nice_my(size(data_PET,1));
        
        data = data_all(ind_sel_data,:);

        opt_rank = app.RankPlotCheckBox.Value;
        hold(app.FigureNeuro, 'on');

        % plot first one only
        for j = 1:size(data_PET,1)
                x = data_PET(j,:);
                y = data(1,:);
                if opt_rank == 1
                    x = tiedrank(x);
                    y = tiedrank(y);
                end
                plot(app.FigureNeuro,x,y,'o','MarkerSize',9,'MarkerFaceColor',colors(j,:),'MarkerEdgeColor','black');
        end
        

        % plot remaining

        for j = 1:size(data_PET,1)
              for i = 2:size(data,1)
                x = data_PET(j,:);
                y = data(i,:);
                if opt_rank == 1
                    x = tiedrank(x);
                    y = tiedrank(y);
                end
                plot(app.FigureNeuro,x,y,'o','MarkerSize',9,'MarkerFaceColor',colors(j,:),'MarkerEdgeColor','black');
              end
        end
        
      for j = 1:size(data_PET,1)
                for i = 1:size(data,1)
                x = data_PET(j,:);
                y = data(i,:);
                if opt_rank == 1
                    x = tiedrank(x);
                    y = tiedrank(y);
                end
                mx = min(x);
                Mx = max(x);
                my = min(y);
                My = max(y);
                limx = max([abs(mx) abs(Mx)]);
                low_x = prctile(x,5);
                high_x = prctile(y,95);
                xfit = mx:0.1:Mx;    
                [p,dev,STATS] = glmfit(x, y);
                [yfit,dlo,dhi] = glmval(p,xfit,'identity',STATS,0.95);
                set(app.FigureNeuro,'fontsize',18);

                DELTA_max = yfit + dlo;
                DELTA_min = yfit - dhi;

                y_area =[DELTA_min', fliplr(DELTA_max')];
                x_area= [xfit, fliplr(xfit)];
                faceAlpha = 0.2;
                plot(app.FigureNeuro, xfit, yfit, 'r-', 'LineWidth', 2,'color',colors(j,:));
                patch(app.FigureNeuro,x_area,y_area,1,'facecolor',colors(j,:),'edgecolor','none','facealpha',faceAlpha);    
                end
        end

        grid(app.FigureNeuro,'on');

        if opt_rank == 1
            xlabel(app.FigureNeuro,'Rank data PET)','FontSize',16);
            ylabel(app.FigureNeuro,'Rank data modality','FontSize',16);
        else
            ylabel(app.FigureNeuro,'Data modality','FontSize',16);
            xlabel(app.FigureNeuro,'Data PET','FontSize',16);
        end
        legend_plot = Rec_list;
        legend(app.FigureNeuro, legend_plot,'AutoUpdate','off','Location','southeast');

    end
   
    
    function [name_save_new] = add_suffix(name_save,options)
    file_part = '';
        switch options(1)
            case 1
                file_part = [file_part '_ESb'];
            case 2
                file_part = [file_part '_ESw'];
            case 3
                file_part = [file_part '_mList1'];
            case 4
                file_part = [file_part '_List1Each'];
            case 5
                file_part = [file_part '_indZ'];
            case 6
                file_part = [file_part '_pwDiff'];
            case 7
                file_part = [file_part '_looList1'];
            case 8
                file_part = [file_part '_List1EachGroupTest'];
        end

        switch options(2)
            case 1
                file_part = [file_part '_Spearman'];
            case 2
                file_part = [file_part '_Pearson'];
            case 3
                file_part = [file_part '_multReg'];
        end

        if options(3) == 1 
            file_part = [file_part '_withExactP'];
        end

        if options(5) == 1
            file_part = [file_part '_withExactSpatialP'];
        end
        time_now = datestr(datetime('now'),'ddmmmyyyy_HHMMSS');
        name_save_new = [name_save file_part '_' time_now];
    end



            

