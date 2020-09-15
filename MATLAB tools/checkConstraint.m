function D_res = checkConstraint(polinomio, val, toll)
% checkConstrain(polinomio, dom, segno, val)
% polinomio: vettore dei coefficienti del polinomio in una variabile
% val: valore presente dopo il segno nel vincolo
% toll: tolleranza assegnata
% D_res: insieme dei domini risultante

% Il dominio della variabile si assume essere in [0,1]
dom(1,1) = 0;
dom(2,1) = 1;

pBernstein = toBernsteinMono(polinomio);
% Porto in forma di Bernstein il polinomio

massimo = max(pBernstein);
minimo = min(pBernstein);
if(massimo > val && minimo > val)
    D_res(1,1) = NaN; % Vincolo mai verificato
    D_res(2,1) = NaN;
elseif(massimo < val+toll && minimo < val+toll)
    D_res(1,1) = dom(1,1); % Vincolo sempre verificato
    D_res(2,1) = dom(2,1);
else
    D = biseziona(polinomio, pBernstein, dom, val, toll);
    
    % Scorro la matrice D per eliminare i NaN che sono stati
    % inseriti in essa dall'algoritmo di bisezione
    indexRes = 1;
    for i = 1 : length(D)
        if(~isnan(D(1,i)))
            D_res(1,indexRes) = D(1,i);
            D_res(2,indexRes) = D(2,i);
            indexRes = indexRes + 1;
        end
    end
end

end

