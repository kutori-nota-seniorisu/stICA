function [ST,locs]=est_spikes_deconv(X,fs,R,nMad,mpd)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:    S - unfused tetanus (velocity)
%           fs - sample rate
%           R - extension factor
%           nMad - multiple of mean absolute deviations
%           mpd - 'MinPeakDistance' in findpeaks
%
% Output:   ST - estimated spike train
%           locs - estimated time instants of the spikes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_iter=100;
Tolx=0.0001;

% Extend
R = R + rem((R+1),2);
[n,N] = size(X);
V = zeros(n*R,N);
for i=1:n
    V((i-1)*R+1:i*R,:)=toeplitz([X(i,1);zeros(R-1,1)],X(i,:));
end
X=V;

% Subtract mean (include in whitening?)
X = bsxfun(@minus, X, mean(X,2));

% Whiten
covarianceMatrix = cov(X');
[E,D] = eig(covarianceMatrix);
whiteningMatrix = sqrt(D)\E';
X = whiteningMatrix*X;

[m,~]=size(X);
w_new=randn(m,1);

n=0;

while true
    % Fixed point algorithm
    w_old=w_new;
    A=mean(tanh(w_old'*X));
    w_new=mean(X.*log(cosh(w_old'*X)),2)-A.*w_old;

    % Normalization
    w_new=w_new/norm(w_new);

    n=n+1;
    if norm(w_new'*w_old-1) < Tolx || n>max_iter
        break;
    end
end

s=w_new'*X;

[~,locs]=findpeaks(s,'MinPeakHeight',nMad*mad(s),'MinPeakDistance',mpd*(fs/1e03));

ST=zeros(size(s));
ST(locs)=1;