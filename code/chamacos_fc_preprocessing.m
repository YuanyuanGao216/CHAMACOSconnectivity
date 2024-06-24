
clc 
clear all
close all

DAP_file = fullfile('..', 'data','chamacos_fNIRS_variables.xlsx');
DAP_table = readtable(DAP_file);
subjlist_from_DAP = cellfun(@str2num,DAP_table{:,1});

%% Load data
datafolder = '/Volumes/Projects/NIRS_Projects/ECHO/y18/derivatives_20220708_1617';
d = dir(datafolder);

d(~[d.isdir])= [];
for ii = numel(d):-1:1
    if strcmp(d(ii).name,'.') || strcmp(d(ii).name,'..')
        d(ii) = [];
    end
end
subjList = {d.name}';
subjList(295:end) = [];    % We don't have DAP data after subj 21963
subjList_number = cellfun(@str2num,subjList(1:end));
for i = 1:length(subjList_number)
    if subjList_number(i) < 1000000
        subjList_number(i) = subjList_number(i)*100+4; 
    end
end
subjList_number(187) = 1509005;

[C,ia] = setdiff(subjList_number,subjlist_from_DAP);
subjList_number(ia) = [];
subjList(ia) = [];
% load the data: 'card'-- card sorting; 'leme'--Sternberg; 'word'--ppt
tasks = {'card';'lmem';'word'}; 

raw_card = cell(numel(subjList),1);
raw_lmem = cell(numel(subjList),1);
raw_word = cell(numel(subjList),1);

for ii = 1:numel(tasks)
     for jj = 1:numel(subjList)
         clear temp data_dir
         data_dir = fullfile(datafolder, subjList{jj}, [subjList{jj} '_' tasks{ii}], [subjList{jj} '_' tasks{ii} '_homer.mat']);
         if isfile(data_dir)
             temp = load(data_dir);
         else
             temp = NaN;
         end
         switch ii
             case 1
                 raw_card{jj} = temp;
             case 2
                 raw_lmem{jj} = temp;
             case 3
                 raw_word{jj} = temp;
         end         
     end        
end

save(fullfile('..', 'data','rawdata.mat'),'raw_card', 'raw_lmem','raw_word', '-v7.3')
save(fullfile('..', 'data','subjList.mat'), 'subjList_number','subjList')
%% ROI

T_proj = readtable(fullfile('..','data','T_proj_clean.xlsx'));
load(fullfile('..','data','Standard_probeInfo.mat'))
x = probeInfo.probes.coords_c2(:,1);
y = probeInfo.probes.coords_c2(:,2);

x_3 = probeInfo.probes.coords_c3(:,1);
y_3 = probeInfo.probes.coords_c3(:,2);
z_3 = probeInfo.probes.coords_c3(:,3);

[roi_labels,IA,IC] = unique(T_proj.roi);
color_array = jet(length(roi_labels));
T_ROI_position = table;
T_ROI_position.roi = roi_labels;
T_ROI_position.roi_number = T_proj.roi_number(IA);
x_mean = zeros(size(roi_labels));
y_mean = zeros(size(roi_labels));
x_mean_3 = zeros(size(roi_labels));
y_mean_3 = zeros(size(roi_labels));
z_mean_3 = zeros(size(roi_labels));
figure

for i = 1:length(roi_labels)
    label = roi_labels{i};
    color = color_array(i,:);
    index_array = find(contains(T_proj.roi, label));

    x_array = x(T_proj{index_array,'ch'});
    y_array = y(T_proj{index_array,'ch'});
    
    plot(x_array, y_array,'LineStyle','none','Marker','o',...
        'MarkerEdgeColor','none','MarkerFaceColor', color, 'MarkerSize',10,'DisplayName',label);
    hold on

    x_mean(i) = mean(x_array);
    y_mean(i) = mean(y_array);
    
    x_array_3 = x_3(T_proj{index_array,'ch'});
    y_array_3 = y_3(T_proj{index_array,'ch'});
    z_array_3 = z_3(T_proj{index_array,'ch'});

    x_mean_3(i) = mean(x_array_3);
    y_mean_3(i) = mean(y_array_3);
    z_mean_3(i) = mean(z_array_3);
end
T_ROI_position.x_mean = x_mean;
T_ROI_position.y_mean = y_mean;
T_ROI_position.x_mean_3 = x_mean_3;
T_ROI_position.y_mean_3 = x_mean_3;
T_ROI_position.z_mean_3 = x_mean_3;
T_ROI_position = sortrows(T_ROI_position, 'roi_number');
writetable(T_ROI_position, fullfile('..','data','T_ROI_position.csv'))


%% ------------------- Data analysis session ---------------- %%
T_proj = readtable(fullfile('..','data','T_proj_clean.xlsx'));
T_ROI_position = readtable(fullfile('..','data','T_ROI_position.csv'));
load(fullfile('..','data','rawdata.mat'));

%% Wisconsin Card sorting: control = 3, match = 2, block_rest = 4, iti_rest = 1
% 1=color, 2=shape, 3=number, 4=control
% matchBlocks = [4,2,4,3,4,2,4,3,4,1,4,2,4,1,4,1,4,3];
fc_card.match = nan(15,15,length(subjList_number));
fc_card.control = nan(15,15,length(subjList_number));
for ii = 1:length(raw_card)
    if ~isstruct(raw_card{ii}) && isnan(raw_card{ii})
        continue
    end
        
    data = raw_card{ii}.hmr;
    evt = data.s;
    hbo = data.hbo;
    fs = data.sampRate;
    badchannel_card = find(data.sciQM.good_combo_link(:,3)<0.7);  % here we use 0.7 as threshold, which is adjustable
    
    % segment data based on evt markers
    evt_blockRest = evt(:,4);
    evt_itiRest = evt(:,1);
    evt_match = evt(:,2);
    evt_control = evt(:,3);
    dur_rest = round(15*fs);
   
    % Match condition
    t = find(evt_match == 1);
    endpoint = [0;find(diff(t)>100)];
    data_match = [];
    for jj = 1:length(endpoint)
            if jj<length(endpoint)
            tp = [t(endpoint(jj)+1):t(endpoint(jj+1))+dur_rest];
            else
                tp = [t(endpoint(jj)+1):t(end)+dur_rest];
            end
        temp_hbo = hbo(tp,:)-mean(hbo(tp,:));% de-centralize the data
        data_match = [data_match;temp_hbo];
    end
    
    
    % Control condition
    t = find(evt_control == 1);
    endpoint = [0;find(diff(t)>100)];
    data_control = [];
    for jj = 1:length(endpoint)
            if jj<length(endpoint)
            tp = [t(endpoint(jj)+1):t(endpoint(jj+1))+dur_rest];
            else
                tp = [t(endpoint(jj)+1):t(end)+dur_rest];
            end
        temp_hbo = hbo(tp,:)-mean(hbo(tp,:));   % de-centralize the data
        data_control = [data_control;temp_hbo];
        clear temp_hbo
    end
    
    
    % ---- Functional connectivity   
    data_match(:,badchannel_card) = nan;
    data_match_roi = zeros(size(data_match,1), height(T_ROI_position));
    data_control(:,badchannel_card) = nan;
    data_control_roi = zeros(size(data_control,1), height(T_ROI_position));
    for j = 1:height(T_ROI_position)
        roi_number = T_ROI_position{j,'roi_number'};
        index = find(T_proj.roi_number == roi_number);
        data_match_roi(:,j) = nanmean(data_match(:,index),2);
        data_control_roi(:,j) = nanmean(data_control(:,index),2);
    end
    fc_match = corr(data_match_roi,data_match_roi);
    fc_control = corr(data_control_roi,data_control_roi);
    
    fc_match(eye(size(fc_match))==1) = nan;
    fc_control(eye(size(fc_control))==1) = nan;
    fc_card.match(:,:,ii) = atanh(fc_match);
    fc_card.control(:,:,ii) = atanh(fc_control);
end
save(fullfile('..','data','fc_card_roi.mat'), 'fc_card')

%% Sternberg task: encoding = 1, rest = 4
load(fullfile('..','data','rawdata.mat'),'raw_lmem');
fc_stern.lmem = nan(15,15,length(raw_lmem));
for ii = 1:length(raw_lmem)
    if ~isstruct(raw_lmem{ii}) && isnan(raw_lmem{ii})
        continue
    end
    data = raw_lmem{ii}.hmr;
    evt = data.s;
    hbo = data.hbo;
    fs = data.sampRate;
    badchannel_lmem = find(data.sciQM.good_combo_link(:,3)<0.7);  % here we use 0.7 as threshold, which is adjustable
    
    % data segmentation
    if size(evt,2) < 4
        fc_stern.lmem(:,:,ii) = nan;
    else
        evt_encode = find(evt(:,1)==1);
        evt_Rest = find(evt(:,4)==1);
        
        if length(evt_encode) < 30
            fc_stern.lmem(:,:,ii) = nan;
        else
            dur_rest = round(15*fs);
            tp = [evt_encode(1):evt_Rest(16)+dur_rest evt_encode(16):evt_Rest(end)+dur_rest];
            data_lmem = hbo(tp,:)-mean(hbo(tp,:));   % de-centralize the data
            % ---- Functional connectivity
            data_lmem(:,badchannel_lmem) = nan;
            data_lmem_roi = zeros(size(data_lmem,1), height(T_ROI_position));
            for j = 1:height(T_ROI_position)
                roi_number = T_ROI_position{j,'roi_number'};
                index = find(T_proj.roi_number == roi_number);
                data_lmem_roi(:,j) = nanmean(data_lmem(:,index),2);
            end
            fc_lmem = corr(data_lmem_roi,data_lmem_roi);
            % exclude bad channels
            fc_lmem(eye(size(fc_lmem))==1) = nan;
            fc_stern.lmem(:,:,ii) = atanh(fc_lmem);
        end
    end
end
save(fullfile('..','data','fc_stern_roi.mat'), 'fc_stern')


%% Pyramids and Palm Trees (PPT) Task: size(control) = 3, word(target test) = 2
load(fullfile('..','data','rawdata.mat'), 'raw_word');
for ii = 1:length(raw_word)
    if ~isstruct(raw_word{ii}) && isnan(raw_word{ii})
        continue
    end
    data = raw_word{ii}.hmr;
    evt = data.s;
    hbo = data.hbo;
    fs = data.sampRate;
    badchannel_word = find(data.sciQM.good_combo_link(:,3)<0.5);  % here we use 0.7 as threshold, which is adjustable
    
    % segment data based on evt markers
    evt_word = evt(:,2);
    evt_size = evt(:,3);
    evt_rest = evt(:,4);
    t_word = find(evt_word == 1);
    t_size = find(evt_size == 1);
    t_rest = find(evt_rest == 1);
    
    if length(t_word) ~= length(t_size)
        fc_ppt.word(:,:,ii) = nan;
        fc_ppt.size(:,:,ii) = nan;

    else
        
        dur_rest = round(10*fs);
        data_word = [];
        for jj = 1:length(t_word)
            if jj<length(t_word)
                tp = t_word(jj):t_rest(2*jj)+dur_rest;
            else
                tp = t_word(jj):t_rest(end)+dur_rest;
            end
            if tp(end) < size(hbo,1)
                temp_hbo = hbo(tp,:)-mean(hbo(tp,:));% de-centralize the data
                data_word = [data_word;temp_hbo];
            end
        end
        % Control condition
        data_size = [];
        for jj = 1:length(t_size)
            tp = t_size(jj):t_rest(2*jj-1)+dur_rest;
            temp_hbo = hbo(tp,:)-mean(hbo(tp,:));   % de-centralize the data
            data_size = [data_size;temp_hbo];
        end
        
        
        % ---- Functional connectivity
        data_word(:,badchannel_word) = nan;
        data_word_roi = zeros(size(data_word,1), height(T_ROI_position));
        data_size(:,badchannel_word) = nan;
        data_size_roi = zeros(size(data_size,1), height(T_ROI_position));
        for j = 1:height(T_ROI_position)
            roi_number = T_ROI_position{j,'roi_number'};
            index = find(T_proj.roi_number == roi_number);
            data_word_roi(:,j) = nanmean(data_word(:,index),2);
            data_size_roi(:,j) = nanmean(data_size(:,index),2);
        end
        fc_word = corr(data_word_roi,data_word_roi);
        fc_size = corr(data_size_roi,data_size_roi);
        % exclude bad channels
        fc_word(eye(size(fc_word))==1) = nan;
        fc_size(eye(size(fc_size))==1) = nan;
        
        %% Put everything into a mat file
        fc_ppt.word(:,:,ii) = atanh(fc_word);
        fc_ppt.size(:,:,ii) = atanh(fc_size);
    end
end
save(fullfile('..','data','fc_ppt_roi.mat'), 'fc_ppt')