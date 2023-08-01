clear all;
close all;
clc;

T = 1;
Ts = 0.0001;

n = [-20:20];
u1 = (n>-T & n<T);
figure();
plot(n,u1);
axis([-7 7 -7 7]);
title('Rff(t)');
set(gcf,'color', 'w');
%2
u2 = (n>-T+2 & n<T+2);
figure();
plot(n,u2);
axis([-20 20 -5 5]);
title('Rff(t-2)');
set(gcf,'color', 'w');
%3
n = -20:Ts:20;
u3 = ((-n/T)-1).* (n> -T & n <= -T/2)+ ((3*n/T)+1).* (n> -T/2 & n <= 0) + (1-(3*n/T)).* (n> 0 & n <= T/2) + ((n/T)-1).* (n> T/2 & n <= T);
figure();
plot(n,u3);
axis([-5 5 -2 2]);
title('Rff(t)');
set(gcf,'color', 'w');
