function [tempi1, tempi2] = confrontoTempi(tentativi)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

load('UnaVariabile')

% Preallocazione per aumentare la velocità
tempi1(tentativi) = 0;
tempi2(tentativi) = 0;

% Prove per l'algoritmo tradizionale
for k = 1 : tentativi
    tstart = tic;
    checkConstraint(pol, val, toll); 
    x = toc(tstart);
    tempi1(k) = x;
    pause(1);
end
n1 = length(tempi1);
media1 = mean(tempi1);
errore1 = std(tempi1)/sqrt(n1);
fprintf('Lo algoritmo classico ha ottenuto: %f +- %f\n', media1, errore1)

% Prove per l'algoritmo di Nataraj
for k = 1 : tentativi
    tstart = tic;
    NatarajCheckConstraint(pol, val, toll);
    x = toc(tstart);
    tempi2(k) = x;
end
n2 = length(tempi2);
media2 = mean(tempi2);
errore2 = std(tempi2)/sqrt(n2);
fprintf('Lo algoritmo di Nataraj ha ottenuto: %f +- %f\n', media2, errore2)

% Grafico degli andamenti risultante
subplot(2,1,1), plot(tempi1), axis([0 1000 0 2])
xlabel('Tentativi')
ylabel('Tempo')
title('Algoritmo tradizionale')
subplot(2,1,2), plot(tempi2), axis([0 1000 0 4*10^-3])
xlabel('Tentativi')
ylabel('Tempo')
title('Algoritmo di Nataraj')


end

