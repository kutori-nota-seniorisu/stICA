function [Data_decomp,Data_array,Data_new,prohibitInd,decompChannelInd] = PreProcess4GUI_v2(Data,parameters,prohibitInd)
%%%%% version 2 was written on 2019-4-13 by CC, which optimized the process order
%%%%% of each step, aiming to clearly input and output the data. 

%%%%% the function was optimized for GUI, which change the input data
%%%%% format from many variables to only one struct variable.

%%%%% first written on 2018-12-19 by CC %%%%%
% INPUT
%   Data: the sEMG data input in a matrix (channels*samples)
%   fsamp: sampling rate
%   electrode: electrode configuration
%   Diff_time: flag determining whether to perform difference in time domain
%   Diff_spatial: flag determining whether to do difference in spatial domain
% OUTPUT
%   Data_new: the filtered data with the same size as Data
%   Data_array: the filtered data with spatial distribution from electrode configuration 


%%%修改：
%24/9/1 从filter修改成filtfilt 
    % 差距不大，使用filter的原因是online反解时加窗下filtfilt边缘效应很强

fsamp = parameters.fsamp;
try
    electrode = parameters.electrodeType;
catch
    electrode = parameters.ElectrodeType;
end
try
    Diff_time = parameters.TimeDifference;
catch
    Diff_time = parameters.timeDifference;
end
try
    Diff_spatial = parameters.spatialDifference;
catch
    Diff_spatial = parameters.SpatialDifference;
end
try
    BPFilter = parameters.bandpassFilter;
catch
    BPFilter = parameters.BandpassFilter;
end
try
    LineFilter = parameters.lineFilter;
catch
    LineFilter = parameters.LineFilter;
end



try
    ChannelFilter = parameters.ChannelFilter;
catch
    ChannelFilter = 0;%?
end
if nargin<3
    prohibitInd = 0;
end

%% re-arrange the Data based on the electrode array configuration
% electrode array set
% 5*13 array, inter-electrode distance 8mm
% back perspective, with adapter below the array
switch electrode
    case 1
        elearr = [4 5 11 10 24 32 34 39 40 49 50 62 61
                  3 6 12 9  23 31 33 38 48 41 51 63 60
                  2 7 13 17 22 30 27 37 47 42 52 64 59
                  1 8 14 18 21 29 26 36 46 43 53 56 58
               NaN 16 15 19 20 28 25 35 45 44 54 55 57];
        elearr = rot90(elearr);
        elearr = rot90(elearr);
        elearr = rot90(elearr);  
        elearr(6,:) = NaN;
    case 2
        % 8*4 array,inter-electrode distance 10mm
        elearr = [8 16 24 32 
                  7 15 23 31 
                  6 14 22 30 
                  5 13 21 29 
                  4 12 20 28 
                  3 11 19 27 
                  2 10 18 26 
                  1 9  17 25  ];
    case 3
        % 8*16 array
        elearr = [8 16 24 32 40 48 56 64
                  7 15 23 31 39 47 55 63
                  6 14 22 30 38 46 54 62
                  5 13 21 29 37 45 53 61
                  4 12 20 28 36 44 52 60
                  3 11 19 27 35 43 51 59
                  2 10 18 26 34 42 50 58
                  1 9  17 25 33 41 49 57];
        elearr = [elearr,elearr+64];
    case 4
        elearr = [4 5 11 10 24 32 34 39 40 49 50 62 61
                  3 6 12 9  23 31 33 38 48 41 51 63 60
                  2 7 13 17 22 30 27 37 47 42 52 64 59
                  1 8 14 18 21 29 26 36 46 43 53 56 58
               NaN 16 15 19 20 28 25 35 45 44 54 55 57];    
    case 5 %左下角1右上角63
        elearr = [8 16 24 32 40 48 56 64
                  7 15 23 31 39 47 55 63
                  6 14 22 30 38 46 54 62
                  5 13 21 29 37 45 53 61
                  4 12 20 28 36 44 52 60
                  3 11 19 27 35 43 51 59
                  2 10 18 26 34 42 50 58
                  1 9  17 25 33 41 49 57];
        % elearr = rot90(elearr);
        % elearr = rot90(elearr);
    case 6
        elearr = [4 5 11 10 24 32 34 39 40 49 50 62 61
                  3 6 12 9  23 31 33 38 48 41 51 63 60
                  2 7 13 17 22 30 27 37 47 42 52 64 59
                  1 8 14 18 21 29 26 36 46 43 53 56 58
               NaN 16 15 19 20 28 25 35 45 44 54 55 57];  
        elearr = elearr'; 
    case 7
        elearr = zeros(11,17);
        for tmpc = 1:17
            elearr(:,tmpc) = (1:11)'+(tmpc-1)*11;
        end
        clear tmpc
    case 8
        elearr = zeros(11);
        for tmpc = 1:11
            elearr(:,tmpc) = (1:11)'+(tmpc-1)*11;
        end
        clear tmpc
    case 9
        elearr = zeros(8);
        for tmpc = 1:8
            elearr(:,tmpc) = [1:8]'+(tmpc-1)*8;
        end
        clear tmpc
    case 10
        elearr = [4 5 11 10 24 32 34 39 40 49 50 62 61
                  3 6 12 9  23 31 33 38 48 41 51 63 60
                  2 7 13 17 22 30 27 37 47 42 52 64 59
                  1 8 14 18 21 29 26 36 46 43 53 56 58
               NaN 16 15 19 20 28 25 35 45 44 54 55 57];
        elearr = rot90(elearr);

    case 11
        elearr = [4 5 11 10 24 32 34 39 40 49 50 62 61
                  3 6 12 9  23 31 33 38 48 41 51 63 60
                  2 7 13 17 22 30 27 37 47 42 52 64 59
                  1 8 14 18 21 29 26 36 46 43 53 56 58
               NaN 16 15 19 20 28 25 35 45 44 54 55 57];  
        elearr = rot90(elearr); 
        elearr = rot90(elearr);
        elearr = rot90(elearr);
    case 12
        % 8*16 array
        elearr = [8 16 24 32 40 48 56 64
                  7 15 23 31 39 47 55 63
                  6 14 22 30 38 46 54 62
                  5 13 21 29 37 45 53 61
                  4 12 20 28 36 44 52 60
                  3 11 19 27 35 43 51 59
                  2 10 18 26 34 42 50 58
                  1 9  17 25 33 41 49 57];
        elearr = [elearr,elearr+64,elearr+128];
    case 13
        % 5*13 array,横着，右上角为缺口
        elearr(1,:) = [12:-1:1,NaN];
        elearr(2,:) = 13:25;
        elearr(3,:) = 38:-1:26;
        elearr(4,:) = 39:51;
        elearr(5,:) = 64:-1:52;
    case 14
        % 8*16 array
        elearr = [8 16 24 32 40 48 56 64
                  7 15 23 31 39 47 55 63
                  6 14 22 30 38 46 54 62
                  5 13 21 29 37 45 53 61
                  4 12 20 28 36 44 52 60
                  3 11 19 27 35 43 51 59
                  2 10 18 26 34 42 50 58
                  1 9  17 25 33 41 49 57];
        elearr = [elearr+128,elearr+64,elearr];
    case 15
        % 8*16 array
        elearr = [7 15 23 31 39 47 55 63
                  6 14 22 30 38 46 54 62
                  5 13 21 29 37 45 53 61
                  4 12 20 28 36 44 52 60
                  3 11 19 27 35 43 51 59
                  2 10 18 26 34 42 50 58
                  1 9  17 25 33 41 49 57];
    case 16
%         elearr = 1:60; % for open-source data
        elearr = reshape([1:60,NaN,NaN,NaN,NaN],8,8);
%     case 17
% %         elearr = 1:60; % for open-source data
%         elearr = reshape([1:60,NaN,NaN,NaN,NaN],8,8);
    case 17
        elearr = zeros(8);
        for tmpc = 1:8
            elearr(:,tmpc) = [1:8]'+(tmpc-1)*8;
        end
        elearr = rot90(elearr); 
        clear tmpc
    case 18%8*8
        elearr = [64,56,48,40,32,24,16,8
                  63,55,47,39,31,23,15,7
                  62,54,46,38,30,22,14,6
                  61,53,45,37,29,21,13,5
                  60,52,44,36,28,20,12,4
                  59,51,43,35,27,19,11,3
                  58,50,42,34,26,18,10,2
                  57,49,41,33,25,17, 9,1];
    case 19
        elearr = [52:64
                  51:-1:39
                  26:38
                  25:-1:13
                  NaN,1:12];
    case 20
        elearr = [1:32];
    case 21
        elearr = 1:16;
    case 22
        elearr = [25:32,NaN,NaN,40:-1:33
                   17:24,NaN,NaN,48:-1:41
                   9:16,NaN,NaN,56:-1:49
                   1:8,NaN,NaN,64:-1:57];
    case 23
        elearr = [1:8:25,NaN,NaN,33:8:57
                  2:8:26,NaN,NaN,34:8:58
                  3:8:27,NaN,NaN,35:8:59
                  4:8:28,NaN,NaN,36:8:60
                  5:8:29,NaN,NaN,37:8:61
                  6:8:30,NaN,NaN,38:8:62
                  7:8:31,NaN,NaN,39:8:63
                  8:8:32,NaN,NaN,40:8:64];
end

%% filtering in time domain
chNum = size(Data,1);
Data_new = Data;
if BPFilter
    LowFreq = 10;   % low frequency of passband 
    HighFreq = 500; % high frequency of passband
    % signalQuality = 5;  % frequency criterion

    % Bandpass filter the vector array Y, pass band LowFreq-HighFreq 
    [Be,Ae] = butter(4,[LowFreq/fsamp*2 HighFreq/fsamp*2]); % EMG band pass filter
    
    for ch = 1:chNum
        Data_new(ch,:) = filtfilt(Be,Ae,Data_new(ch,:));%24/9/1 从filter修改成filtfilt
    end
end

% remove abnormal frequency components
if LineFilter
    % use the function from Holobar
%         Data_new = chooseChannels(Data_new,fsamp,1,0,50);

%     % use the bandstop filter
%     LowFreq = 49;
%     HighFreq = 51;
%     [Be,Ae] = butter(4,[LowFreq/fsamp*2 HighFreq/fsamp*2],'stop'); % EMG band pass filter

    % Optimized the line filter: change the bandstop filter at 50 Hz to the
    % comb filter -- by CC at 2019-10-21
    fo = 50; % power frequency
    q = 20;
    bw = (fo/(fsamp/2))/q;
    [Be,Ae] = iircomb(round(fsamp/fo),bw,'notch');%%%梳状滤波
    
    for ch = 1:chNum
        Data_new(ch,:) = filtfilt(Be,Ae,Data_new(ch,:));%24/9/1 从filter修改成filtfilt
    end

    %%%%%
    % 使用Movi时可能有111Hz，需要111Hz滤波
    % fo = 111; % power frequency
    % q = 10;
    % bw = (fo/(fsamp/2))/q;
    % [Be,Ae] = iircomb(round(fsamp/fo),bw,'notch');%%%梳状滤波
    % 
    % for ch = 1:chNum
    %     Data_new(ch,:) = filter(Be,Ae,Data_new(ch,:));
    % end
    %%%%%
end

% perform difference in time domain
if Diff_time
    Data_new = diff(Data_new')';
end
        
%% channel selection depending on the amplitude feature
if ChannelFilter
    if prohibitInd==0   
        AmpCrit = zeros(1,chNum);
        for ch = 1:chNum
            AmpCrit(ch) = rms(Data(ch,:));
        end
        tmpave = mean(AmpCrit(find(AmpCrit~=0)));
        tmpstd = std(AmpCrit(find(AmpCrit~=0)));
        AmpCrit(find(AmpCrit > tmpave+3*tmpstd)) = 0;
        AmpCrit(find(AmpCrit <= max(tmpave-3*tmpstd,0))) = 0;
        Ind = find(AmpCrit==0);
        Data_new(Ind,:) = 0; 
        prohibitInd = Ind;
    else        
        Data_new(prohibitInd,:) = 0;
    end
end
    
%% filtering in spatial domain
Data_array = cell(size(elearr));
for r = 1:size(elearr,1)
    for c = 1:size(elearr,2)
        ch = elearr(r,c);        
        if ~isnan(ch) && ~ismember(ch,prohibitInd)
%             sum(abs(Data_new(ch,:)))~=0 
            Data_array{r,c} = Data_new(ch,:);
        end
    end
end

if Diff_spatial
    % Data_array = bipolarfilter(Data_array);
    Data_array{1,13}=Data_array{1,12};
    for i=1:5
        for j=1:12
            % tmpmuap{1,1}=zeros(1,50);
            Data_tmp{i,j}=Data_array{i,j+1}-Data_array{i,j};
        end
    end
    Data_array=Data_tmp;
end

Data_decomp = [];
decompChannelInd = [];
for c = 1:size(Data_array,2)
    for r = 1:size(Data_array,1)
        if ~isempty(Data_array{r,c})
            Data_decomp(end+1,:) = Data_array{r,c};
            decompChannelInd(end+1) = elearr(r,c);
        end
    end
end

