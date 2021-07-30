function varargout = JuSpace(varargin)
% JuSpace MATLAB code for JuSpace.fig
%      JuSpace, by itself, creates a new JuSpace or raises the existing
%      singleton*.
%
%      H = JuSpace returns the handle to a new JuSpace or the handle to
%      the existing singleton*.
%
%      JuSpace('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JuSpace.M with the given input arguments.
%
%      JuSpace('Property','Value',...) creates a new JuSpace or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JuSpace_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JuSpace_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JuSpace

% Last Modified by GUIDE v2.5 05-Jul-2021 15:31:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JuSpace_OpeningFcn, ...
                   'gui_OutputFcn',  @JuSpace_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end


if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before JuSpace is made visible.
function JuSpace_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JuSpace (see VARARGIN)

% Choose default command line output for JuSpace
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
dir_tool= fileparts(which('JuSpace'));
% dir_spm = which('spm');
% [path,dir] = fileparts(dir_spm);
dir_atlas = fullfile(dir_tool,'atlas','m_labels_Neuromorphometrics.nii');
set(handles.atlaslist,'String', dir_atlas);

dir_PET = fullfile(dir_tool,'PETatlas');
global files_PET
files_PET = select_con_maps_forfMRI_my(dir_PET,'PETatlas','.*.nii');

for i = 1:length(files_PET)
    [path,file] = fileparts(files_PET{i});
    files{i} = file;
end

set(handles.PETlist,'String',files);
% UIWAIT makes JuSpace wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = JuSpace_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in compute.
function compute_Callback(hObject, eventdata, handles)
% hObject    handle to compute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dir_save = get(handles.dirsave,'String');

% get options
opt_comp = find([handles.Opt_comp.Children.Value]==1);
% opt_comp = 1 --> es between
% opt_comp = 2 --> es within
% opt_comp = 3 --> mean list 1
% opt_comp = 4 --> list 1 each
% opt_comp = 5 --> ind z-score list 1 to list 2
% opt_comp = 6 --> pair-wise difference list 1 to list 2
% opt_comp = 7 --> ind loo z-scores from list 1
% opt_comp = 8 --> list 1 each compares against null distribution of
% correlation coefficients
opt_ana = find([handles.Ana_opt.Children.Value]==1);
opt_perm = get(handles.opt_perm,'Value');
opt_perm_spatial = get(handles.opt_perm_spatial,'Value');
opt_auto = get(handles.opt_autocorr_T1,'Value');
Nperm = str2num(get(handles.nperm,'String'));
options = [opt_comp; opt_ana; opt_perm; opt_auto; opt_perm_spatial];

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

% get files
list1 = get(handles.filelist1,'String');
list2 = get(handles.filelist2,'String');

if isempty(list1)
    f = msgbox('Please select files in list 1 for correlation with PET');
    return
end

if and(~ismember(opt_comp,[3,4,7,8]),isempty(list2))
    f = msgbox('This option requires files to be selected for list 2');
    return
end

if isempty(dir_save)
    f = msgbox('Please select first a directory to save the results');
    return
end

if or(opt_comp == 2,opt_comp== 6)
    if length(list1)~=length(list2)
    f = msgbox('Pair-wise design requires same number of images in list 1 and 2');
    return
    end
end

if and(opt_comp == 5,length(list2)<3)
    f = msgbox('This option requires more than 2 files in list 2');
    return
end

if and(opt_comp == 2,length(list1)<3)
    f = msgbox('This option requires more than 2 files in list 1 and 2');
    return
end

if and(opt_comp == 7,length(list1)<3)
    f = msgbox('This option requires more than 2 files in list 1');
    return
end


global files_PET
aa = get(handles.PETlist,'Value');
str_PET = get(handles.PETlist,'String');
str_PET = str_PET(aa);
for i = 1:length(str_PET);
   tt = regexp(str_PET{i},'_','Split');
   Rec_list{i}=tt{1};
end
filesPET = files_PET(aa);
atlas = char(get(handles.atlaslist,'String'));

name_save = get(handles.namesave,'String');
time_now = datestr(datetime('now'),'ddmmmyyyy_HHMMSS');

file_save = fullfile(dir_save,[name_save file_part '_' time_now '.mat']);

if options(1)<4
    image_save = fullfile(dir_save,[name_save file_part '_' time_now '.nii']);
    [res,p_all,stats,data,D1,D2, data_PET,Resh,T1] = compute_DomainGauges(list1,list2,filesPET,atlas, options,image_save);
else
    [res,p_all,stats,data,D1,D2, data_PET,Resh,T1] = compute_DomainGauges(list1,list2,filesPET,atlas, options);
end
%
opt_for_perm = [1,2,5,6];
opt_for_spat_perm = [3, 4, 8];

if options(3)==1 && ismember(options(1),opt_for_perm)% && options(2)~=3
    disp('Computing exact p-value');
   [p_exact,dist_r] = compute_exact_pvalue(D1,D2,data_PET,res,Nperm,options,T1);
   Resh(:,end+1) = [{'p_exact'}; num2cell_my(p_exact')];
end

if options(5)==1 && ismember(options(1),opt_for_spat_perm)
    disp('Computing exact spatial p-value')
    [p_exact,dist_r] = compute_exact_spatial_pvalue(D1,data_PET,atlas,res,Nperm,options,filesPET, T1);
    Resh(:,end+1) = [{'p_exact'}; num2cell_my(p_exact')];
end


switch options(1)
    case 1
        ff = 'Effect Size';
    case 2
        ff = 'Effect Size';
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
    case 7
        ff = 'Z-score loo';
    case 8
        ff = 'List 1 all against null';
end


if exist('dist_r','var')
    save(file_save,'res','p_all','stats','data','Resh','D1','D2','data_PET','list1','list2','filesPET','atlas','options','dist_r','p_exact','Nperm');
else
    save(file_save,'res','p_all','stats','data','Resh','D1','D2','data_PET','list1','list2','filesPET','atlas','options');
end

res = stats.res_ind;
ind_plot =[];
for i = 1:size(res,1)
    for j = 1:size(res,2)
       ind_plot(end+1,1) = j;
    end
end

vals_plot= res; %cell2num_my(Resh(2:end,3));
for i = 1:length(filesPET)
   all_i = vals_plot(:,i);%vals_plot(ind_plot==i);
   
   size_y_xi = sum(ind_plot==i);
   y_dist = abs(vals_plot(ind_plot==i) - mean(vals_plot(ind_plot==i)));
   y_dist_inv = abs(y_dist -max(y_dist))./max(abs(y_dist -max(y_dist)));
   if isnan(y_dist_inv)
       y_dist_inv = rand(size(y_dist_inv));
   end
   
   x_n_n(ind_plot==i) = ind_plot(ind_plot==i) + 1.*(rand(size_y_xi,1)-0.5).*y_dist_inv;
   all_ii = all_i(~isinf(all_i));
   m_x(i) = mean(all_ii);
   std_x(i,1) = 1.95996.*std(all_ii)./sqrt(length(all_ii));
%    if opt_comp>=4
%        
%    end
end
% plot specs
marker_color = [0 0.7 1];
bar_color = generate_colors_nice_my(length(filesPET));
dot_edge = [0 0.4 1];
error_color = [0.2 0.3 0.4];
error_thickness = 3;
marker_size = 20;
%_______

%-----------------------
% generate bar plots
%-----------------------
% hold('off');
set(handles.res_summary,'Data',Resh);
if length(filesPET)<3
    h1 = figure('Position',[400 400 500 500]);
else
    h1 = figure('Position',[400 400 length(filesPET).*120 500]);
end

h2 = bar(1:length(filesPET),diag(m_x),0.9,'stacked');
%colormap(bar_color)
% bar(handles.plot_bars, 1:length(filesPET),m_x,0.9,'EdgeColor',bar_edge);
hold('on');

ax = gca;
set(ax,'XTickLabel',Rec_list,'FontSize',16);
vv = vals_plot';
vv2 = vv(:);

for i = 1:length(filesPET)
    set(h2(i),'facecolor',bar_color(i,:));
    if opt_comp>=4
        plot(x_n_n(ind_plot==i),vv2(ind_plot==i),'o','MarkerSize',8,'MarkerFaceColor',bar_color(i,:),'MarkerEdgeColor',error_color,'LineWidth',1);
    end
end
if opt_comp>=4
        h2 = errorbar(1:length(filesPET),m_x,std_x,'.');
        title('Error bars represent 95% confidence interval','FontSize',12);
        set(h2,'LineWidth',error_thickness,'Color',error_color,'MarkerSize',marker_size,'MarkerEdgeColor',error_color);
end



% if opt_ana<3
%    ylim([-1 1]);
% end


% t = findall(gcf,'Type','Axes');
% pp = FDR(p_all,0.05);
% for i = 1:length(p_all)
%     if p_all(i)<=pp
%         if opt_comp<4
%             text(t(end),i,max(m_x)+0.2,'*','FontSize',20);
%     end
% end

switch opt_ana
    case 1
        ylabel(ax,'Fisher''s z (Spearman rho)','fontweight','bold');
    case 2
        ylabel(ax,'Fisher''s z (Pearson r)','fontweight','bold');
    case 3
        ylabel(ax,'Beta coeffiecient','fontweight','bold');
end
savefig(h1,fullfile(dir_save,['Bar_' name_save file_part '_' time_now '.fig']));
print(h1,fullfile(dir_save,['Bar_' name_save file_part '_' time_now '.png']),'-dpng','-r300');
%-----------------------
% end generate bar plots
%-----------------------

%-----------------------
% generate scatter plots (only for less than 10 lines)
%-----------------------
% if [size(data,1).*size(data_PET,1)]<10

if opt_ana<3
%     n_color = 1./(n_draw+1);
    n_color = 1./(size(data,1)+1);
    n_rows = ceil(size(data_PET,1)./3);
    hh = figure('Position',[400 0 900 n_rows.*300]);
    subplot(n_rows,3,1);
    n_draw = 31;
    for j = 1:size(data_PET,1)
        nn = 1;
        colors_all = [];
        color_i = [1 0.5 0];
%         legend_all = {};
        subplot(n_rows,3,j);
%         if opt_ana == 1
%             xlabel(['Relative rank ' Rec_list{j}],'FontSize',12);
%         else
            xlabel(['Data ' Rec_list{j}],'FontSize',12);
%         end
        hold('on');
            for i = 1:size(data,1)
                    nn = nn + 1;
                    if nn < n_draw+1
                        x = data_PET(j,:);
                        y = data(i,:)';
%                         if opt_ana == 1
%                             [temp,x]  = ismember(x,unique(x));
%                             [temp,y]  = ismember(y,unique(y));
%                         end
                        plot(x,y,'o','MarkerSize',6,'MarkerFaceColor',color_i,'MarkerEdgeColor','black');
                        colors_all = [colors_all; color_i];
                        color_i(1) = color_i(1)-n_color;
                        color_i(3) = color_i(3)+n_color;
%                         legend_all{end+1} = [ff num2str(i) ' with ' Rec_list{j}];
                    end

                end

%         if size(data,1)==1
%             legend(legend_all,'AutoUpdate','off','Location','southoutside');
%         else
%         end
        h1 = lsline;
        set(h1,'LineWidth',3);
        for i = 1:length(h1)
           h1(i).Color = colors_all(i,:); 
        end

    end
        h1=subplot(n_rows,3,1);
        h2=subplot(n_rows,3,3);
        h3=subplot(n_rows,3,n_rows.*3 - 2);
%         h4=subplot(n_rows,3,n_rows.*3);
        p1=get(h1,'position');
        p2=get(h2,'position');
        p3=get(h3,'position');
%         p4=get(h4,'position');
        height=p3(2) + (p1(2)-p3(2)) + p1(4)./2;
        width=p3(1)+ (p2(1)-p3(1)+ p3(3)./2);
        h5=axes('position',[p3(1) p3(2).*0.7 width height],'visible','off'); 
        h5.XLabel.Visible='on';
        h5.YLabel.Visible='on';
        axes(h5);
%         if opt_ana == 1
%             ylabel('Relative rank modality','FontSize',20);
%         else
            ylabel('Data Modality','FontSize',20);
%         end
        xlabel('PET','FontSize',20);
    
    
    
%     suplabel('Data PET','x','FontSize',20);
%     suplabel('Data Modality','y','FontSize',20);
    savefig(hh,fullfile(dir_save,['Line_' name_save file_part '_' time_now '.fig']));
    print(hh,fullfile(dir_save,['Line_' name_save file_part '_' time_now '.png']),'-dpng','-r300');
end
%-----------------------
% end generate scatter plots
%-----------------------

% --- Executes on button press in selectfiles1.
function selectfiles1_Callback(hObject, eventdata, handles)
% hObject    handle to selectfiles1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
files1 = cellstr(spm_select(Inf,'image','Select files for spatial correlation'));
if ~isemptycell(files1)
    set(handles.filelist1,'String',files1);
end


% --- Executes on button press in selectfiles2.
function selectfiles2_Callback(hObject, eventdata, handles)
% hObject    handle to selectfiles2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
files2 = cellstr(spm_select(Inf,'image','Select files for spatial correlation'));
if ~isemptycell(files2)
    set(handles.filelist2,'String',files2);
end


% --- Executes on button press in selectatlas.
function selectatlas_Callback(hObject, eventdata, handles)
% hObject    handle to selectatlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dir_cd = cd;
path = fileparts(which('JuSpace'));
path_atlas = fullfile(path,'atlas');
cd(path_atlas);
files_atlas = cellstr(spm_select(1,'image','Select atlas file you would like to use'));
cd(dir_cd);

if ~isemptycell(files_atlas)
    set(handles.atlaslist,'String',files_atlas);
end


% --- Executes during object creation, after setting all properties.
function filelist1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filelist1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function atlaslist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to atlaslist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function filelist2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filelist2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in selectdirectory.
function selectdirectory_Callback(hObject, eventdata, handles)
% hObject    handle to selectdirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dir_save = uigetdir(cd, 'Select a directory to save the results');
if ~(dir_save == 0)
    set(handles.dirsave,'String',dir_save);
end


% --- Executes during object creation, after setting all properties.
function dirsave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dirsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function namesave_CreateFcn(hObject, eventdata, handles)
% hObject    handle to namesave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function PETlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PETlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
clear handles;
clear all


% --- Executes on selection change in PETlist.
function PETlist_Callback(hObject, eventdata, handles)
% hObject    handle to PETlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns PETlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PETlist


% --- Executes on selection change in filelist1.
function filelist1_Callback(hObject, eventdata, handles)
% hObject    handle to filelist1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filelist1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filelist1


% --- Executes on selection change in filelist2.
function filelist2_Callback(hObject, eventdata, handles)
% hObject    handle to filelist2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filelist2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filelist2



function namesave_Callback(hObject, eventdata, handles)
% hObject    handle to namesave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of namesave as text
%        str2double(get(hObject,'String')) returns contents of namesave as a double



function dirsave_Callback(hObject, eventdata, handles)
% hObject    handle to dirsave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dirsave as text
%        str2double(get(hObject,'String')) returns contents of dirsave as a double


% --- Executes on button press in opt_regr_PET.
function opt_regr_PET_Callback(hObject, eventdata, handles)
% hObject    handle to opt_regr_PET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of opt_regr_PET


% --- Executes on button press in opt_corr_PET.
function opt_corr_PET_Callback(hObject, eventdata, handles)
% hObject    handle to opt_corr_PET (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of opt_corr_PET


% --- Executes on button press in opt_perm.
function opt_perm_Callback(hObject, eventdata, handles)
% hObject    handle to opt_perm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of opt_perm



function nperm_Callback(hObject, eventdata, handles)
% hObject    handle to nperm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nperm as text
%        str2double(get(hObject,'String')) returns contents of nperm as a double


% --- Executes during object creation, after setting all properties.
function nperm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nperm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function Opt_comp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Opt_comp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in opt_autocorr_T1.
function opt_autocorr_T1_Callback(hObject, eventdata, handles)
% hObject    handle to opt_autocorr_T1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of opt_autocorr_T1


% --- Executes on button press in opt_perm_spatial.
function opt_perm_spatial_Callback(hObject, eventdata, handles)
% hObject    handle to opt_perm_spatial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of opt_perm_spatial
