function [toll, tempi1, tempi2] = confrontoTempi3()
%CONFRONTOTEMPI3 Summary of this function goes here
%   Detailed explanation goes here

load('UnaVariabile')
format long
%toll = [0.1 0.09 0.08 0.07 0.06 0.05 0.04 0.03 0.02 0.01 10^-3 10^-4 10^-5];
toll = [10^-1 10^-2 10^-3 10^-4 10^-5 10^-6 10^-7 10^-8 10^-9 10^-10];

% Prove per l'algoritmo tradizionale
for k = 1 : length(toll)
    tstart = tic;
    checkConstraint(pol, val, toll(k));
    x = toc(tstart);
    tempi1(k) = x;
    pause(1);
end

% Prove per l'algoritmo di Nataraj
for k = 1 : length(toll)
    tstart = tic;
    NatarajCheckConstraint(pol, val, toll(k)); 
    x = toc(tstart);
    tempi2(k) = x;
    pause(1);
end

% Grafico degli andamenti risultante
subplot(2,1,1), semilogx(toll, tempi1)%, axis([0 0.1 0 1])
xlabel('Tolleranza')
ylabel('Tempo')
title('Algoritmo tradizionale')
subplot(2,1,2), semilogx(toll, tempi2)%, axis([0 0.1 0 1])
xlabel('Tolleranza')
ylabel('Tempo')
title('Algoritmo di Nataraj')


end


