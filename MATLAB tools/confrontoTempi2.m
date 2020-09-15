function [tempi1, tempi2] = confrontoTempi2(tentativi)

load('UnaVariabile')

% Preallocazione per aumentare la velocità
tempi1(tentativi) = 0;
tempi2(tentativi) = 0;



% Prove per l'algoritmo tradizionale
for k = 1 : tentativi
    tic; checkConstraint(pol, val, toll); tempi1(k) = toc;
end
n1 = length(tempi1);
media1(1) = mean(tempi1);
errore1(1) = std(tempi1)/sqrt(n1);
fprintf('Lo algoritmo classico ha ottenuto: %f +- %f\n', media1, errore1)

% Prove per l'algoritmo di Nataraj
for k = 1 : tentativi
    tic; NatarajCheckConstraint(pol, val, toll); tempi2(k) = toc;
end
n2 = length(tempi2);
media2(1) = mean(tempi2);
errore2(1) = std(tempi2)/sqrt(n2);
fprintf('Lo algoritmo di Nataraj ha ottenuto: %f +- %f\n', media2, errore2)



pol = [0 1 -10 10];
val = 10;
% Prove per l'algoritmo tradizionale
for k = 1 : tentativi
    tic; checkConstraint(pol, val, toll); tempi1(k) = toc;
end
n1 = length(tempi1);
media1(2) = mean(tempi1);
errore1(2) = std(tempi1)/sqrt(n1);
fprintf('Lo algoritmo classico ha ottenuto: %f +- %f\n', media1, errore1)

% Prove per l'algoritmo di Nataraj
for k = 1 : tentativi
    tic; NatarajCheckConstraint(pol, val, toll); tempi2(k) = toc;
end
n2 = length(tempi2);
media2(2) = mean(tempi2);
errore2(2) = std(tempi2)/sqrt(n2);
fprintf('Lo algoritmo di Nataraj ha ottenuto: %f +- %f\n', media2, errore2)




% Grafico degli andamenti risultante
subplot(2,1,1), plot(tempi1)
xlabel('Tentativi')
ylabel('Tempo')
title('Algoritmo tradizionale')
subplot(2,1,2), plot(tempi2), axis([0 1000 0 0.008])
xlabel('Tentativi')
ylabel('Tempo')
title('Algoritmo di Nataraj')


end

