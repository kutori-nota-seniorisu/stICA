% 尝试绘制六参数模型曲线。实现方式已修改无误
clc; clear; close all;

% ms
Thc = 30;
Tc = 60;
Thr = 125;
Ttw = 240;
Tlead = 0;

Dmax = 11;

% 240ms
t = 0:Ttw;

e = exp(1);
c1 = log(2)*Tc / (Thc - Tc + Tc*log(Tc/Thc));
c2 = log(2)*Tc / (Thr - Tc + Tc*log(Tc/Thr));

P1 = (t - Tlead) / Tc;
P2 = 1 + exp(2*e*(P1-1));
P3 = (t - 0.5*(Ttw+Thr)) / (Ttw-Thr);

f1 = Dmax * (P1.^c1 .* exp(c1-c1*P1) + (P2-1) .* P1.^c2 .* exp(c2-c2*P1));
f2 = P2 + P2.*exp(4*e*(P3));

F = f1 ./ f2;

figure;
plot(F);

