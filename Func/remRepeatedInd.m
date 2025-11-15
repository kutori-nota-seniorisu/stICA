function NewPulseInd = remRepeatedInd(PulseSeq,PulseInd,limit)
% NewPulseInd = remRepeatedInd(PulseSeq,PulseInd,limit);
%
% Checks the interpulse distance of pulses contained in "PulseInd".
% If the interpulse distances drops below the "limit" the powers of 
% two investigated subsequent pulses are compared and the pulse 
% with the smallest power is left out.  
% 
% Inputs: 
%   PulseSeq - tested sequence of pulses
%   PulseInd - time positions of tested pulses.
%   limit - the lower limit of allowed interpulse distance
% Outputs:
%   NewPulseInd - new (filtered) time positions of pulses


r = length(PulseInd);
PulseInd = sort(PulseInd);
tmp = zeros(1,length(PulseInd));
NewPulseInd = tmp;
pulse_step = 1;
NewPulseInd(1) = -inf;

if (r<2)
    NewPulseInd = PulseInd;
    return 
end 

k = 1;
while k < r
    while (k<r && (PulseInd(k)-NewPulseInd(pulse_step))<limit)    % find the first position
        k = k+1;
    end 
    first = k;
    tmp_step = 1;
    tmp(tmp_step) = PulseSeq(PulseInd(k));
    while (k<r && (PulseInd(k+1)-PulseInd(first))<limit)
        % compare the next index(first+1) value and the present(first)
        % index value, if the difference is over the limit, no operation
        % will be executed. if the difference drops below the limit, the
        % next index(first+1) will be compared with first index. the
        % loop stops until finding a index distance over the limit
        tmp_step = tmp_step+1;
        tmp(tmp_step) = PulseSeq(PulseInd(k+1));
        k = k+1;
    end 
    % compare the index and find the maximum value (highest power)
    [~,i] = max(tmp(1:tmp_step));
    pulse_step = pulse_step+1;
    NewPulseInd(pulse_step) = PulseInd(first-1+i);
    k = k+1;
end 

if (PulseInd(end)-PulseInd(end-1)>limit)
    pulse_step = pulse_step+1;
    NewPulseInd(pulse_step) = PulseInd(end);
end 
NewPulseInd = NewPulseInd(2:pulse_step);  % the index is the array number

