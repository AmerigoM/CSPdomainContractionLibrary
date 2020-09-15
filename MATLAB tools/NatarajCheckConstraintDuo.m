function M_res = NatarajCheckConstraintDuo(polinomio, dom1, dom2, val, toll)
% Restrizione domini di variabili per vincoli bidimensionali
% polinomio: matrice dei coefficienti del polinomio in due variabili
% dom1: matrice 2x1 del dominio della prima variabile del polinomio
% dom2: matrice 2x1 del dominio della seconda variabile del polinomio
% val: valore che costituisce il vincolo
% toll: tolleranza assegnata per l'accettazione dell'intervallo
% M_res: matrice 4xN dove ogni colonna rappresenta un "pezzo" del dominio
%        della superficie tridimensionale. I primi due valori sono per la
%        prima variabile mentre gli altri due sono per la seconda variabile

% EURISTICA: taglio il dominio della variabile con dominio maggiore

dim = size(polinomio);
% Ricavo le dimensioni della matrice (deve essere "minima")

dom = cat(2,dom1,dom2);
polB = toBernsteinDuo(polinomio, dom);
% Porto il polinomio in due variabili in base di Bernstein

massimo = max(max(polB));
minimo = min(min(polB));
% Calcolo il massimo e il minimo coefficiente di Bernstein

if(massimo > val+toll && minimo > val+toll)
% Vincolo mai verificato
    M_res(1,1) = NaN;
    M_res(2,1) = NaN;
    M_res(3,1) = NaN;
    M_res(4,1) = NaN;
elseif(massimo < val+toll && minimo < val+toll)
% Vincolo sempre verificato
    M_res = dom1;
    M_res = [M_res; dom2];
else
% Non è possibile una deduzione immediata
    D_res = NatarajBisezionaDuo(polB, dim(1)-1, dim(2)-1, dom1, dom2, val, toll);
    
    % Elimino tutti i NaN inseriti all'interno di D_res
    i = 1;
    dimention = size(D_res);
    for k = 1 : dimention(2)
        if(~isnan(D_res(1,k)))
            M_res(:,i) = D_res(:,k);
            i = i+1;
        end
    end
    
end


end

