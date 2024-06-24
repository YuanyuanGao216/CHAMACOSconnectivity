% this the original code I got from Eric
clc 
clear all
close all
%% Hardcoded block order
% 1=color, 2=shape, 3=number, 4=control
matchBlocks = [4,2,4,3,4,2,4,3,4,1,4,2,4,1,4,1,4,3];



%% Load data
datafolder = '/Volumes/EricSSD/CHAMACOS/derivatives_20220708_1617';
d = dir(datafolder);

d(~[d.isdir])= [];
for ii = numel(d):-1:1
    if strcmp(d(ii).name,'.') || strcmp(d(ii).name,'..')
        d(ii) = [];
    end
end
subjList = {d.name}';
subjList(455) = [];    % last one is not a subj folder, delete it

% load the data: 'card'-- card sorting; 'leme'--Sternberg; 'word'--ppt
tasks = {'card';'lmem';'word'}; 

raw_card = [];
raw_lmem = [];
raw_word = [];

for ii = 1:numel(tasks)
     missing_sub{ii} = {};
    
     for jj = 1:numel(subjList)
         clear temp data_dir
         data_dir = [datafolder filesep subjList{jj} filesep [subjList{jj} '_' tasks{ii}] filesep [subjList{jj} '_' tasks{ii} '_homer.mat']];
         if isfile(data_dir)
             temp = load(data_dir);
             switch ii
                 case 1
                     raw_card = [raw_card;temp];
                 case 2
                     raw_lmem = [raw_lmem;temp];
                 case 3
                     raw_word = [raw_word;temp];
             end
         else
             missing_sub{ii} = [missing_sub{ii};subjList{jj}];
             warning([subjList{jj} '_' tasks{ii} ' do not exsit, please check!'])
         end
         
         
     end        
end

return


%% ------------------- Data analysis session ---------------- %%
load('rawdata.mat');

%% Wisconsin Card sorting: control = 3, match = 2, block_rest = 4, iti_rest = 1
% 1=color, 2=shape, 3=number, 4=control
%matchBlocks = [4,2,4,3,4,2,4,3,4,1,4,2,4,1,4,1,4,3];
channels = 36;

clear hbo_match hbo_control fc_card badchannel_card
for ii = 1:numel(raw_card)
    
    clear data_match data_control fc_match fc_control
    data = raw_card(ii).hmr;
    evt = data.s;
    hbo = data.hbo;
    fs = data.sampRate;
    badchannel_card{ii} = find(data.sciQM.good_combo_link(:,3)<0.7);  % here we use 0.7 as threshold, which is adjustable
    
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
        clear temp_hbo
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
    
    
    %% ---- Functional connectivity    
    fc_match = corr(data_match,data_match);
    fc_control = corr(data_control,data_control);
    % exclude bad channels
    fc_match(badchannel_card{ii},:) = 0;
    fc_match(:,badchannel_card{ii}) = 0;
    fc_match = fc_match.*~eye(size(fc_match));

    fc_control(badchannel_card{ii},:) = 0;
    fc_control(:,badchannel_card{ii}) = 0;
    fc_control = fc_control.*~eye(size(fc_control));

    %% Put everything into a mat file
    hb_card.hbo_match{ii} = data_match;
    hb_card.hbo_control{ii} = data_control;
    fc_card.match(:,:,ii) = atanh(fc_match);
    fc_card.control(:,:,ii) = atanh(fc_control);
      
    
end

for ii = 1:numel(badchannel_card)
   ct(ii) = length(badchannel_card{ii}); 
    
end
length(find(ct>=20))

sum(ct)/36/numel(badchannel_card)

hist(ct)



%% Sternberg task: encoding = 1, rest = 4

clear hbo_lmem fc_stern badchannel_lmem
bad_evt_lmem = [];
for ii = 1:numel(raw_lmem)
    
    clear data_lmem temp_hbo fc_lmem
    data = raw_lmem(ii).hmr;
    evt = data.s;
    hbo = data.hbo;
    fs = data.sampRate;
    badchannel_lmem{ii} = find(data.sciQM.good_combo_link(:,3)<0.7);  % here we use 0.7 as threshold, which is adjustable
    
    %% data segmentation
    if size(evt,2) < 4
        hb_lmem.hbo_lmem{ii} = [];
        fc_stern.lmem(:,:,ii) = 0;
        bad_evt_lmem = [bad_evt_lmem;ii];
    else
        evt_encode = find(evt(:,1)==1);
        evt_Rest = find(evt(:,4)==1);
        
        if length(evt_encode) < 30
            hb_lmem.hbo_lmem{ii} = [];
            fc_stern.lmem(:,:,ii) = 0;
            bad_evt_lmem = [bad_evt_lmem;ii];
        else
            dur_rest = round(15*fs);
            data_lmem = [];
            tp = [evt_encode(1):evt_Rest(16)+dur_rest evt_encode(16):evt_Rest(end)+dur_rest];
            
            temp_hbo = hbo(tp,:)-mean(hbo(tp,:));   % de-centralize the data
            data_lmem = [data_lmem;temp_hbo];
            %% ---- Functional connectivity
            fc_lmem = corr(data_lmem,data_lmem);
            % exclude bad channels
            fc_lmem(badchannel_lmem{ii},:) = 0;
            fc_lmem(:,badchannel_lmem{ii}) = 0;
            fc_lmem = fc_lmem.*~eye(size(fc_lmem));
            
            %% Put everything into a mat file
            hb_lmem.hbo_lmem{ii} = data_lmem;
            fc_stern.lmem(:,:,ii) = atanh(fc_lmem);
        end
    end
end

for ii = 1:numel(badchannel_lmem)
   ct(ii) = length(badchannel_lmem{ii}); 
    
end
length(find(ct>=30))

sum(ct)/36/numel(badchannel_lmem)



%% Pyramids and Palm Trees (PPT) Task: size(control) = 3, word(target test) = 2


clear hb_ppt fc_ppt hbo_word hbo_size badchannel_word
bad_evt_word = [];
for ii = 1:numel(raw_word) 
    
    clear data_word data_size fc_word fc_size hbo
    data = raw_word(ii).hmr;
    evt = data.s;
    hbo = data.hbo;
    fs = data.sampRate;
    badchannel_word{ii} = find(data.sciQM.good_combo_link(:,3)<0.5);  % here we use 0.7 as threshold, which is adjustable
    
    % segment data based on evt markers
    evt_word = evt(:,2);
    evt_size = evt(:,3);
    evt_rest = evt(:,4);
    t_word = find(evt_word == 1);
    t_size = find(evt_size == 1);
    t_rest = find(evt_rest == 1);
    
    if length(t_word) ~= length(t_size)
        hb_ppt.hbo_word{ii} = [];
        hb_ppt.hbo_size{ii} = [];
        fc_ppt.word(:,:,ii) = 0;
        fc_ppt.size(:,:,ii) = 0;
        bad_evt_word = [bad_evt_word;ii];

    else
        % Word condition
        dur_rest = round(10*fs);
        data_word = [];
        for jj = 1:length(t_word)
            
            if jj<length(t_word)
                tp = [t_word(jj):t_rest(2*jj)+dur_rest];
            else
                tp = [t_word(jj):t_rest(end)+dur_rest];
                
            end
            if tp(end) < size(hbo,1)
                temp_hbo = hbo(tp,:)-mean(hbo(tp,:));% de-centralize the data
                data_word = [data_word;temp_hbo];
            end
            clear temp_hbo
        end
        
        
        % Control condition
        data_size = [];
        for jj = 1:length(t_size)
            
            tp = [t_size(jj):t_rest(2*jj-1)+dur_rest];
            
            temp_hbo = hbo(tp,:)-mean(hbo(tp,:));   % de-centralize the data
            data_size = [data_size;temp_hbo];
            clear temp_hbo
        end
        
        
        %% ---- Functional connectivity
        fc_word = corr(data_word,data_word);
        fc_size = corr(data_size,data_size);
        % exclude bad channels
        fc_word(badchannel_word{ii},:) = 0;
        fc_word(:,badchannel_word{ii}) = 0;
        fc_word = fc_word.*~eye(size(fc_word));
        
        fc_size(badchannel_word{ii},:) = 0;
        fc_size(:,badchannel_word{ii}) = 0;
        fc_size = fc_size.*~eye(size(fc_size));
        
        %% Put everything into a mat file
        hb_ppt.hbo_word{ii} = data_word;
        hb_ppt.hbo_size{ii} = data_size;
        fc_ppt.word(:,:,ii) = atanh(fc_word);
        fc_ppt.size(:,:,ii) = atanh(fc_size);
    end
    
end
for ii = 1:numel(badchannel_word)
   ct(ii) = length(badchannel_word{ii}); 
    
end
length(find(ct>=30))

sum(ct)/36/numel(badchannel_word)


%% ----------------- Compute graph theory measure -------------------- %%
roi = {[1 3 5];[2 6];[4 7]; [10 11]; [8 12]; [9 13]; [14 15 17]; [16 18]; ...
       [19 21 24]; [23 25]; [20 22 26]; [27 28]; [29 33 34]; [30 31 35]; [32 36]} ;
roi_num = numel(roi);


for ii = 1:numel(badchannel_card)
   ct(ii) = length(badchannel_card{ii}); 
    
end
sum(ct)/36/numel(badchannel_word)

hist(ct)


group = cell(2*length(nonzeros(ge_card_match)),1);
group(1:length(nonzeros(ge_card_match))) = {'match'};
group(length(nonzeros(ge_card_match))+1:end) = {'control'};
%%
up = 1;
low = 0.1;
step = 0.05;

%% WSC card sorting
[~,ge_card_match] = chomsFC_RL(fc_card.match, badchannel_card, up, low, step);
[~,ge_card_control] = chomsFC_RL(fc_card.control, badchannel_card, up,low, step);

[h,p]=ttest(nonzeros(ge_card_match),nonzeros(ge_card_control))
[p,h]=signrank(nonzeros(ge_card_match),nonzeros(ge_card_control))

figure,
h1 = histogram(nonzeros(ge_card_match)), hold on
h2 = histogram(nonzeros(ge_card_control))
%h1.Normalization = 'probability';
h1.BinWidth = 0.01;
%h2.Normalization = 'probability';
h2.BinWidth = 0.01;
xlabel('global efficiency')
ylabel('count number')
legend('card sorting match','card sorting control')

figure(1),
vs = violinplot([nonzeros(ge_card_match) nonzeros(ge_card_control)]);
ylabel('global efficiency');
xlim([0 3]);

%% Sternberg (lmem)
[~,ge_stern] = chomsFC_RL(fc_stern.lmem, badchannel_lmem, up,low, step);


%% PPT word
[~,ge_ppt_word] = chomsFC_RL(fc_ppt.word, badchannel_word, up, low, step);
[~,ge_ppt_size] = chomsFC_RL(fc_ppt.size, badchannel_word, up, low, step);
[h,p]=ttest(nonzeros(ge_ppt_word),nonzeros(ge_ppt_size))
[p,h]=signrank(nonzeros(ge_ppt_word),nonzeros(ge_ppt_size))



figure,
h1 = histogram(nonzeros(ge_stern))
h1.BinWidth = 0.01;
xlabel('global efficiency')
ylabel('count number')
legend('Sternberg encoding')

figure,
h1 = histogram(nonzeros(ge_ppt_word)), hold on
h2 = histogram(nonzeros(ge_ppt_size))
%h1.Normalization = 'probability';
h1.BinWidth = 0.01;
%h2.Normalization = 'probability';
h2.BinWidth = 0.01;
xlabel('global efficiency')
ylabel('count number')
legend('PPT semantic word','PPT size/control')







