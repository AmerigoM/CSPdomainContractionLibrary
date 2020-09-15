function res = toBernsteinMono(C)
% Function per la conversione di un polinomio da Forma di Potenze
% a Forma di Bernstein: res = toBernsteinMono(C)
% C: vettore dei coefficienti del polinomio in una variabile
%    (dove il termine più a sinistra è il coefficiente di x^0)
% res: vettore dei coefficienti del polinomio in base di Bernstein

degree = length(C) - 1;
% Calcolo il grado del polinomio che gli è stato passato

for k = 1 : length(C)
% Calcolo i primi coefficienti con la formula ricavata dallo sviluppo di
% Taylor
    B(1,k) = C(k)/nchoosek(degree, k-1);
end

index = 1;
% Indice per calcolare i restanti elementi

for k = 2 : length(C)
% Calcolo le restanti soluzioni (anche parziali) inserendole in una matrice
    for j = 1 : length(C) - index
        B(k,j) = B(k-1,j) + B(k-1,j+1);
    end
    index = index + 1;
end

for k = 1 : length(B)
% Scrivo gli effettivi coefficienti di Bernstein nella variabile da
% ritornare "res"
    res(k) = B(k, 1);
end

end

