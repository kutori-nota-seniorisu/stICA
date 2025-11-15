function [Cost,newPulses,LowLimit] = evaluateIPT(PulseJitt,Pulse,Fsample,LowLimit,NoEpochs,mode)
% function [Cost,newPulses]=evaluateIPT(PulseJitt,Pulse,Fsample,LowLimit,NoEpochs,mode)
%
% Finds and evaluates the optimal MU firing pattern
% INPUTS:
%   PulseJitt - tested pulse sequence
%   Pulse - preselected pulses
%   Fsample - sampling frequency
%   LowLimit - relative limit of allowed impulse height w.r.t. maximal pulse height (possible values lie on the interval [0,1]);
%   NoEpochs - number of consecutive epochs of PulseJitt pulse sequence (the epochs length is calculated automatically)
%   mode - temporal parameter for debugging purposes
%
% OUTPUTS:
%   Cost - costs of identified MU firing pattern,heuristically estimated scores of reconstructed MU discharge patterns (the lower the score,the better the reconstruction)
%   newPulses - identified MU firing pulses;

if nargin > 5
    plotMode = 1;
else 
    plotMode = 0;
end 

% the frequency band of discharge rate -- for experimental signals -
% 2021-11-3
FreqMin = 2;
FreqMax = 50;
PulseDistance = round(Fsample/100);

%%% usually for simulated signals -- 2021-11-3
% FreqMin = 1;
% FreqMax = 35;
% PulseDistance = 25;
% PulseDistance = 50;

step = 0.05;

outliersTsh = 2;
if isempty(PulseJitt)
    if isempty(Pulse)
        Cost = inf;
        newPulses = [];
        return ;
    else 
        MaxNoPulses = FreqMax*(max(1,max(Pulse))-max(1,min(Pulse)))/Fsample;
    end 
else 
    MaxNoPulses = FreqMax*length(PulseJitt)/Fsample;
    % calculate the maximum pulse number in a whole pulse train 
end 
% CostPerPulse = 1/(MaxNoPulses+eps);

Cost = inf;
newPulses = Pulse;
% tmpCostFun = 0;
maxPulse = max(PulseJitt);
tmpLowLimit = LowLimit;

if plotMode > 0
    hfig = figure;set(hfig,'Position',[1,1,1350,970]);
end 

% find a optimized k to achieve a lowest cost value 
for k = 0.9:-step:tmpLowLimit
    if isempty(PulseJitt)
        tmpInd = Pulse;
        if plotMode > 0
            subplot(2,1,1);plot(tmpInd,ones(size(tmpInd)),'or');
        end         
    else
        tmpInd = find(PulseJitt > k*maxPulse);
        tmpInd = remRepeatedInd(PulseJitt,tmpInd,PulseDistance);
        if plotMode > 0       
            figure(hfig);hold off;
            subplot(2,1,1);
            plot(PulseJitt);hold on;
            plot([1,length(PulseJitt)],[k*maxPulse,k*maxPulse],'r');
            plot(tmpInd,PulseJitt(tmpInd),'or');
        end
    end
 
%     % evaluate way from CC
%     if (length(tmpInd) > 5)
%         tmpCostFun = CoVISI({tmpInd},Fsample,[0,max(tmpInd)/Fsample]);
%         if Cost > tmpCostFun
%             Cost = tmpCostFun;
%             newPulses = tmpInd;
%             LowLimit = k;
%         end
%     end
    % origin evaluate way from Dario
    if (length(tmpInd) > 5)
        tmp2 = diff(tmpInd);
        p = polyfit(1:length(tmp2),tmp2,0);   % polynomial curve fitting (0-roder, actually p is a constant equalling the mean of tmp2)
        iIPI = polyval(p,1:length(tmp2));  % calculate the polynomial value 
%         figure;
%         plot(1:length(tmp2),tmp2,'*');
%         hold on;plot(1:length(tmp2),iIPI);
        
        % cut outliers (just for missed firings) 
        tmp = find(abs(tmp2-iIPI)./abs(iIPI)>outliersTsh);
        tmpInd2 = setdiff(1:length(tmp2),tmp);
        if length(tmpInd2) > 6
            tmp2 = tmp2(tmpInd2(3:end-2));
        else
            tmp2 = tmp2(tmpInd2);
        end
        % end cut outliers
        %%
        if (length(tmp2) > 2)
            p = polyfit(1:length(tmp2),tmp2,2);
            iIPI = polyval(p,1:length(tmp2));
            iFreq = Fsample ./ iIPI; % instaneous discharge frequency of each discharge

            Ffir_mean = median(iFreq);
            Ffir_min = min(iFreq);
            Ffir_max = max(iFreq);
            Ffir_res = median(abs(tmp2-iIPI) ./ abs(iIPI));

            if length(tmp2) < 6
                tmpCostFun = 1;
            else
                tmpCostFun = 0;
            end
            
            if Ffir_mean < FreqMin || Ffir_mean > FreqMax
                tmpCostFun = tmpCostFun+100;
            end
            if Ffir_min < FreqMin || Ffir_min > FreqMax
                tmpCostFun = tmpCostFun+100;
            end
            if Ffir_max < FreqMin || Ffir_max > FreqMax
                tmpCostFun = tmpCostFun+100;
            end
            if isempty(PulseJitt)
                MaxNoPulses = Ffir_mean*(max(1,max(Pulse))-max(1,min(Pulse)))/Fsample;
            else
                MaxNoPulses = Ffir_mean*(max(1,max(tmpInd))-max(1,min(tmpInd)))/Fsample;
            end
            CostPerPulse = 1/MaxNoPulses;

%             NoLargeIPI = length(find(Fsample ./ tmp2 < 4));
%             NoSmallIPI = length(find(Fsample ./ tmp2 > 35));

            tmpCostFun = (tmpCostFun+4.0*Ffir_res+1*abs(1-length(tmpInd)*CostPerPulse))/5;

            if plotMode > 0
                title(['Cost = ',num2str(tmpCostFun),'   1-length(tmp2)*CostPerPulse = ',num2str(abs(1-length(tmpInd)*CostPerPulse))]);
                hold off;
                subplot(2,1,2);
                plot(Fsample ./ tmp2,'b*','LineWidth',3,'MarkerSize',7);hold on;plot(Fsample ./ iIPI,'r');
                hold on;
                title(['Ffir_{min} = ',num2str(Ffir_min),'  Ffir_{mean} = ',num2str(Ffir_mean),'   Ffir_{max} = ',num2str(Ffir_max),'  Ffir_{res} = ',num2str(Ffir_res) ...
                    ,'     MaxNoPulses = ',num2str(MaxNoPulses),'  No. Identified pulses = ',num2str(length(tmpInd))]);
%                 pause;
            end
            
            if Cost > tmpCostFun
                Cost = tmpCostFun;
                newPulses = tmpInd;
                LowLimit = k;
%             elseif 3*Cost < tmpCostFun
            end
        end 
    end
end

% if plotMode > 0
%     figure(hfig);close(hfig);
% end 

