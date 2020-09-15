function res = DeCasteljau(pBernstein, point)
% Implementazione dell'algoritmo di De Casteljau ad una dimensione
% pBersntein: vettore dei coefficienti di Bernstein a partire
%             dal coefficiente di B_0^n
% point: punto in cui avviene la valutazione
% res: risultato della valutazione

% Calcolo la dimensione della matrice delle soluzioni parziali
lung = length(pBernstein);

% Inizializzo la matrice delle soluzioni parziali
M = zeros(lung,lung);

% Passo base
for k = 1 : lung
    M(1,k) = pBernstein(k);
end

% Indice per calcolo piramidale nella matrice
index = 1;

% Passo induttivo
for i = 2 : lung
    for j = 1 : lung - index
        M(i,j) = M(i-1,j)*(1-point) + M(i-1,j+1)*point;
    end
    index = index + 1;
end

% Estrazione del risultato finale
res = M(lung,1);

end

