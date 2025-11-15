function DecompResults = IPTExtraction_gCKC4GUI_v3(Sig,parameters)
%%%%%% optimized by CC %%%%%
% calculating the CorrSig without the first and last n points (n=exFactor)
%%%%%% version 2
% extract spikes with the k-means clustering
%%%%%% optimized for EMGDecoder at 2019-3-7 %%%%%%
% the main optimization:
% 		○ 把输入和输出变量改成了结构体
% 		○ 优化了输出变量，方便后续在线解码
%%%%%% first written at 2019-2-21 %%%%%
% this version derives from the v1
% the main optimization:
% 		○ 融合了elimination的步骤
% 		○ 将MUST match的sensitivity阈值改成了0.5，且对含有较少脉冲的MUST匹配进行了优化
% 		○ 优化了pnr的计算代码（只有在IPT被留存的时候才会进行计算）
%       ○ 去掉了w迭代结束后两次extraction的过程
tic;
fsamp = parameters.fsamp;
exFactor = parameters.extendingFactor;
IterationNo_MU = parameters.iterationNumMU;
CostFcn = parameters.costFcn;
drawMode = 0;
iterationNo_w = parameters.iterationNumW;
iterateNoMod = 20;%之前是20
SIL_threshold = 0.8;%博士论文中是0.9，一般0.8
% spikeExtractionFlag = parameters.spikeExtractionFlag;

MUPulses = {}; % the motor unit spike trains
SILs = []; % the confidence indictor of MUSTs
IPTs = []; % the signals from which the MUST is extracted
PNRs = []; % the pulse-to-noise ratio of MUSTs
W = []; % the separation vectors
Centroids = []; % the cluster centroids of each MUST
debugpara.output = []; % the flag indicating the time when each MUST is extracted
debugpara.n0_raw = [];
debugpara.n0 = [];

ActIndSdt = 3; % the threshold of activation level 异常值一般去除3sigma
% iterationNo_w = 45;
stepsize = 1;
sAlpha = 0.1;%学习率
% iterateNoMod = 20;
%% extending
% while median(abs(Sig(:)))<10
%     Sig = 10*Sig;
% end
% Sig(find(sum(Sig,2)==0),:) = [];
eSig = extend(Sig,exFactor);
eSig = eSig(:,exFactor+1:end-exFactor);
%% correlation matrix
% CorrSig = eSig(:,exFactor+1:end-exFactor)*eSig(:,exFactor+1:end-exFactor)'/length(eSig(:,exFactor+1:end-exFactor));
CorrSig = eSig*eSig'/length(eSig);%??没有减去均值（白化）
invCorrSig = pinv(CorrSig);%%%Cxx-1
debugpara.invCorrSig = invCorrSig;
debugpara.eSig = eSig;
%% activation index
ActIndWhole = calcActivityIndex(eSig,invCorrSig,5000);%得到1*length的激活度
% ActIndResid = ActIndWhole;
ActIndResid = fftshift(fft(ActIndWhole));%fft后将零频率分量移动到中心，之后只留下频率分量为正负10%以内的信号
% extract frequency information lower than 10% fsamp (keep the frequency compenent < 204.8 Hz)
ActIndResid([1:round(prctile(1:length(ActIndResid),40)),round(prctile(1:length(ActIndResid),60)):end]) = 0;
ActIndResid = real(ifft(ifftshift(ActIndResid)));

% remove the first and last n points (may be affected by extending)%去掉前后的被扩展影响的点数
edgepoint = exFactor*10;
ActIndResid([1:edgepoint,end-edgepoint+1:end]) = 0;

% remove the abnormal amplitude value
ActIndResid(find(ActIndResid<mean(ActIndWhole)-ActIndSdt*std(ActIndWhole))) = 0;%此处加不加edgepoint区别不大，找异常小值
tmpInd = find(ActIndResid>mean(ActIndWhole(edgepoint+1:end-edgepoint))+ActIndSdt*std(ActIndWhole(edgepoint+1:end-edgepoint)));
prohibitedInd = tmpInd;
while ~isempty(tmpInd)
    tmpInd = tmpInd([1,find(diff(tmpInd)>1)+1]);%移除连续的索引
    for p = 1:length(tmpInd)
        [minInd,~,maxInd] = findHillLim(ActIndResid,tmpInd(p));%找到局部峰值的最大最小索引，之后清除(??为什么不和100行一样)
        ActIndResid(minInd:maxInd) = 0;
    end
    tmpInd = find(ActIndResid>mean(ActIndWhole)+ActIndSdt*std(ActIndWhole));%再次寻找异常大点
    prohibitedInd = [prohibitedInd,tmpInd];
end
prohibitedInd = unique(prohibitedInd);%移除
prohibitedInd = sort(prohibitedInd);%从小到大排序，得到去除的点数号

siglen = size(eSig,2);

if drawMode
    figure1 = figure('color',[1,1,1]);
    set(figure1,'DoubleBuffer','on');
end

fprintf('    Spike extraction ...\n');
%% iteration for IPT extraction
for run_IPT = 1:IterationNo_MU
    % the initial step
    [~,n0] = max(ActIndResid);
    debugpara.n0_raw(run_IPT) = n0;

    [minInd,~,maxInd] = findHillLim(ActIndResid,n0);
    ActIndResid(minInd:maxInd) = 0;
    prohibitedInd = [prohibitedInd, [minInd:maxInd]];

    c_tx = mean(eSig(:,n0),2);
    c_tx = c_tx/sqrt(sum(c_tx.^2)); % normalization

    t_new = c_tx'*invCorrSig*eSig;%%%MUAPT估计,\hat(S_j),c_tx'*invCorrSig是分离向量的转置
    t_new([1:edgepoint,end-edgepoint+1:end]) = 0;
    t_old = zeros(size(t_new));

    similarityFlag = 0;
    %% iteration for separation vector
    for run_c = 1:iterationNo_w
        similarity = sum(abs(t_old-t_new))/sum(abs(t_new));%初始为1
        % break the loop depending on the iteration similarity condition(两次similarity小于0.01，break，即两次估计的误差相差不大)
        if similarity<0.01
            similarityFlag = similarityFlag+1;
            if similarityFlag>=2
                break;
            end
        end
        t_old = t_new;

        % the patial derivative of the cost function
        a1 = 1;
        if CostFcn == 1
            weights = log(1+t_new.^2);
        elseif CostFcn == 2
            weights = log(1+t_new.^4);
        elseif CostFcn == 3
            weights = t_new.^2;
        elseif CostFcn == 4
            weights = 2*a1^2*t_new./(1+a1^2*t_new.^2);
        elseif CostFcn == 5
            weights = abs(t_new);
        else
            weights = t_new.*log(1+a1^2*t_new.^2)-2*t_new+2/a1*atan(a1*t_new);
        end

        if ~mod(run_c,iterateNoMod)
            for run_dc = 1:5
                t_new([1:edgepoint+1,end-edgepoint+1:end]) = 0;
                t_new = t_new/max(t_new);
                tT = abs(t_new).*t_new;
                [~,compInd,~] = evaluateIPT(tT(edgepoint+1:end-edgepoint),[],fsamp,0.0,1);
                c_tx = sum(eSig(:,compInd+edgepoint),2);
                c_tx = c_tx/realsqrt(sum(realpow(c_tx,2)));
                t_new = (c_tx'*invCorrSig)*eSig;
            end
            tT = abs(t_new).*t_new;
            if  -min(t_new) > max(t_new)
                tT(find(t_new>0)) = 0;
            else
                tT(find(t_new<0)) = 0;
            end
            tT = abs(tT);
            %             figure,plot(tT);
            %             hold on;
            %             plot(compInd,tT(compInd),'x');

            [compInd,spikeSIL,C] = spikeExtraction(tT(edgepoint+1:end-edgepoint)',fsamp);
            compInd = compInd+edgepoint;

            %             pnr = 10*log10(mean(tT(compInd).^2)/mean(tT(setdiff([1:length(tT)],compInd)).^2)); % CC for pulse-noise ratio

            if spikeSIL>SIL_threshold
                % remove the these pulse points from the index
                for p = 1:length(compInd)
                    [minInd,~,maxInd] = findHillLim(ActIndResid,compInd(p));
                    ActIndResid( minInd:maxInd ) = 0;
                end
                MUPulses{end+1} = compInd;
                PNRs(end+1) = 10*log10(mean(tT(compInd).^2)/mean(tT(setdiff([1:length(tT)],compInd)).^2)); % CC for pulse-noise ratio
                SILs(end+1) = spikeSIL;
                IPTs(end+1,:) = tT;
                Centroids(end+1,:) = C;

                W(:,end+1) = c_tx; % the separation vector of this IPT
                debugpara.output(:,end+1) = 2;
                debugpara.n0(end+1) = n0;
            end
            clear p dS dTP num1 num2 LowLimit tmpInd minInd maxInd tT

            tmpWeights = zeros(size(weights));
            tmpWeights(compInd) = weights(compInd);
            weights = tmpWeights;
            clear tmpWeights
        end

        weights(prohibitedInd) = 0;
        % the increment of c_tx in the iteration loop
        dc = zeros(size(eSig,1),1);
        for k = edgepoint+1:siglen-edgepoint
            dc = dc+weights(k)*eSig(:,k);
        end
        dc = dc/sqrt(sum(dc.^2));%ctx的增量
        % the increment of t in the iteration loop
        dt = dc'*invCorrSig*eSig;
        dt([1:edgepoint,end-edgepoint+1:end]) = 0;

        % the cost function
        if CostFcn == 1
            oldJ = mean((log(1+a1^2*t_new.^2)-2).*t_new+2/a1.*atan(a1*t_new)); % oldJ derives from t, intergrate of weights, the cost function
        elseif CostFcn == 2
            oldJ = mean(t_new.*log(1+t_new.^4)-4*t_new+1/2*2^(1/2)*log((t_new.^2+2^(1/2)*t_new+1)./...
                (t_new.^2-2^(1/2).*t_new+1))+2^(1/2)*atan(2^(1/2)*t_new+1)+2^(1/2)*atan(2^(1/2)*t_new-1));
        elseif CostFcn == 3
            oldJ = mean(1/3*t_new.^3);
        elseif CostFcn == 4
            oldJ = mean(log(1+(a1*t_new).^2));
        elseif CostFcn == 5
            oldJ = mean((log(1+a1^2*t_new.^2)-2).*t_new+2/a1.*atan(a1*t_new));
        else
            oldJ = mean((1/2*log(1+a1^2*t_new.^2)-3/2).* t_new.^2+2/a1*t_new.*atan(a1*t_new)-1/2/...
                a1^2*log(1+a1^2*t_new.^2)-1/2/a1^2);
        end
        mode = sign(oldJ);
        %更新ctx和tnew
        c_tx = c_tx+mode*stepsize*sAlpha*dc;%??加还是减有说法吗
        clen = sqrt(sum(c_tx.^2));
        c_tx = c_tx/clen;
        t_new = (t_new+mode*stepsize*sAlpha*dt)/clen;%？？同217
        %     t_new([1:edgepoint,end-edgepoint+1:end]) = 0;
        %         t_test = c_tx'*invCorrSig*eSig;

        if drawMode
            tT = t_new.^2;
            if  -min(t_new) > max(t_new)
                tT(find(t_new>0)) = 0;
            else
                tT(find(t_new<0)) = 0;
            end
            [compInd,spikeSIL] = spikeExtraction(tT(edgepoint+1:end-edgepoint)',fsamp);
            compInd = compInd+edgepoint;

            figure(figure1);
            subplot(2,1,1);
            plot([1:length(tT)]/fsamp,tT,compInd/fsamp,tT(compInd),'o');
            %             hold on;
            %             plot(compInd,tT(compInd),'o');
            xlabel('Time(s)','FontSize',14);
            title(['Decomp run: ',num2str(run_IPT),' - SIL: ',num2str(round(spikeSIL*100)/100)],'FontSize',16);
            axis tight;

            subplot(2,1,2);
            plot(round(length(tT)/2)/fsamp+[1:fsamp/2]/fsamp,tT(round(end/2+[1:fsamp/2])));
            xlabel('Time(s)','FontSize',14);
            title(['Iteration k = ' num2str(run_c) ', t^{(k+1)}-t^{(k)}=',num2str(similarity)],'FontSize',14);
            axis tight;
            drawnow;
            saveas(gcf,['C:\Users\kk200\Desktop\draw\MU' num2str(run_IPT) 'Iteration' num2str(run_c)],'jpg');
            pause(0.1);

        end
    end
    %% IPT extraction
    tT = t_new.^2;
    if  -min(t_new) > max(t_new)
        tT(find(t_new>0)) = 0;
    else
        tT(find(t_new<0)) = 0;
    end
    %     [CostFunComp,compInd,LowLimit] = evaluateIPT(tT(edgepoint+1:end-edgepoint),[],fsamp,0);
    %     compInd = compInd+edgepoint;

    [compInd,spikeSIL,C] = spikeExtraction(tT(edgepoint+1:end-edgepoint)',fsamp);

    %     tmptT = tT(edgepoint+1:end-edgepoint);
    %     noiseInd = setdiff([1:length(tmptT)],compInd);
    %     C1 = mean(tmptT(compInd));
    %     C2 = mean(tmptT(noiseInd));
    %
    %     D = [];
    %     D(:,1) = abs(tmptT-C1);
    %     D(:,2) = abs(tmptT-C2);
    %
    %     flag = ones(length(tmptT),1);
    %     flag(compInd) = -1;
    %     SILs = mean(((D(:,1)-D(:,2)).*flag)./max(D,[],2));

    compInd = compInd+edgepoint;

    [~,n02] = max(-t_new);
    [minInd,~,maxInd] = findHillLim(ActIndResid,n02);
    prohibitedInd = unique([prohibitedInd,minInd:maxInd]);

    tT2 = tT/max(tT);
    if spikeSIL>SIL_threshold
        % remove the these pulse points from the index
        for p = 1:length(compInd)
            [minInd,~,maxInd] = findHillLim(ActIndResid,compInd(p));
            ActIndResid( minInd:maxInd ) = 0;
        end

        MUPulses{end+1} = compInd;
        PNRs(end+1) = 10*log10(mean(tT(compInd).^2)/mean(tT(setdiff([1:length(tT)],compInd)).^2)); % CC for pulse-noise ratio
        SILs(end+1) = spikeSIL;
        IPTs(end+1,:) = tT;
        W(:,end+1) = c_tx; % the separation vector of this IPT
        Centroids(end+1,:) = C;
        debugpara.output(:,end+1) = 2;
        debugpara.n0(end+1) = n0;
    end
    clear p dS dTP num1 num2 LowLimit tmpInd minInd maxInd tT
    if drawMode
        figure(figure1);
        if spikeSIL>SIL_threshold
            hold off;plot( [1:length(tT2)]/fsamp,tT2);
            hold on;plot(compInd/fsamp,tT2(compInd),'go');hold off;
            axis tight;xlabel('time (s)','FontSize',14);
            title(['Decomp. run: ',num2str(IterationNo_MU),' - SIL: ',num2str(round(spikeSIL*100)/100)],'FontSize',16);
            drawnow;
        else
            hold off;plot( [1:length(tT2)]/fsamp,tT2);
            hold on;plot(compInd/fsamp,tT2(compInd),'ro');hold off;
            axis tight;xlabel('time (s)','FontSize',14);
            title(['Decomp. run: ',num2str(run_c),' - SIL: ',num2str(round(spikeSIL*100)/100)],'FontSize',16);
            drawnow;
        end
        saveas(gcf,['C:\Users\kk200\Desktop\pics\MU' num2str(run_IPT) 'Iteration' num2str(run_c)],'jpg');
        pause;
    end
end


if drawMode
    try
        close(figure1);
    catch
    end
end
% fprintf([num2str(toc) ' s \n']);
% tic;
fprintf('    Removing duplications: ');
[~,Ind] = sort(SILs,'descend');
k = 0;
newMUPulses = {};
newInd = [];
dIPI = round(fsamp*0.0005);
dataSource = 'exp';
switch dataSource
    case 'sim'
        lim = round(fsamp*0.1);
        fprintf('configuration of duplication elimination: sim. ...\n');
    case 'exp'
        lim = round(fsamp*0.1);
        fprintf('configuration of duplication elimination: exp. ...\n');
end

% matchThresh = 0.5;
while ~isempty(Ind)
    k = k+1;
    cand = Ind(1);

    for i = 2:length(Ind)
        if length(MUPulses{Ind(1)})/length(MUPulses{Ind(i)})<0.3 || length(MUPulses{Ind(i)})/length(MUPulses{Ind(1)})>3.3
            continue;
        end
        r = RoA(MUPulses{Ind(1)},MUPulses{Ind(i)},lim,dIPI);
        if r>0.3
            cand(end+1) = Ind(i);
        end
    end
    tmpSIL = SILs(cand);
    [~,i] = max(tmpSIL);
    newMUPulses{k} = MUPulses{cand(i)};
    newInd(k) = cand(i);
    [~,tmpInd] = setdiff(Ind,cand);
    tmpInd = sort(tmpInd);
    Ind = Ind(tmpInd);
end

% [~,newInd2] = remRepMUST(MUPulses,SILs,dIPI,lim);

newSILs = SILs(newInd);
newPNRs = PNRs(newInd);
newIPTs = IPTs(newInd,:);
newCtx = W(:,newInd);
newCentroids = Centroids(newInd,:);

debugpara.n0 = debugpara.n0(newInd);

DecompResults.MUPulses = newMUPulses;
DecompResults.SILs = newSILs;
DecompResults.IPTs = newIPTs;
DecompResults.PNRs = newPNRs;
DecompResults.Ctx = newCtx;
DecompResults.Centroids = newCentroids;
DecompResults.ActIndWhole = ActIndWhole;
DecompResults.ActIndResid = ActIndResid;

DecompResults.debugpara = debugpara;
fprintf([num2str(toc) ' s \n']);
try
    DecompResults.W = (DecompResults.Ctx'*invCorrSig)';
catch
    DecompResults.W = [];
end