function [polB1, polB2] = NatarajRicalcolaBernstein(degree, pBernstein)
% Function che ricalcola i coefficienti di Bernstein
% degree: grado del polinomio
% pBernstein: coefficienti di Bernstein del polinomio, da sfruttare per
%             effettuare il ricalcolo


% Inserisco i primi coefficienti (noti) nella matrice
for i = 1 : degree + 1
    M(1,i) = pBernstein(i);
    % Matrice delle soluzioni parziali
end

for k = 2 : degree + 1
    for i = 1 : degree + 1
        if(k > i) % Formula di Nataraj
            M(k,i) = M(k-1,i);
        else % k <= i
            M(k,i) = (0.5*M(k-1,i-1)) + (0.5*M(k-1,i));
        end
    end
end


for i = 1 : degree + 1
    % Calcolo gli effettivi coefficienti per la prima metà dell'intervallo
    polB1(i) = M(degree+1,i);
    % Calcolo i coefficienti per la seconda metà dell'intervallo per
    % simmetria
    polB2(i) = M(degree+1-(i-1),degree+1);
end


end
