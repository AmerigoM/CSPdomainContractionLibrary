function D_res = biseziona(polinomio, pBernstein, dom, val, toll)
% Function che applica l'algoritmo di bisezione per restringere il dominio
% polinomio: vettore dei coefficienti del polinomio in una variabile
% pBernstein: vettore dei coefficienti di polinomio in base di Bernstein
% dom: vettore 2x1 del dominio della variabile del polinomio
% val: valore presente dopo il segno nel vincolo
% toll: tolleranza assegnata
% D_res: insieme dei domini risultante

pol = polinomio;
massimo = max(pBernstein);
minimo  = min(pBernstein);

if(massimo > val && minimo > val)
    D_res(1,1) = NaN; % Vincolo mai verificato
    D_res(2,1) = NaN;
elseif(massimo < val+toll && minimo < val+toll)
    D_res(1,1) = dom(1,1); % Vincolo sempre verificato
    D_res(2,1) = dom(2,1);
else
    middle = (dom(1,1) + dom(2,1))/2;
    % Calcolo il punto medio
    
    dom1(1,1) = dom(1,1);
    dom1(2,1) = middle;
    [~, polB1] = ricalcolaBernstein(pol, dom1(1,1), dom1(2,1));
    D1 = biseziona(pol, polB1, dom1, val, toll);
    % Ricorro sulla prima metà dell'intervallo
    
    dom2(1,1) = middle;
    dom2(2,1) = dom(2,1);
    [~, polB2] = ricalcolaBernstein(pol, dom2(1,1), dom2(2,1));
    D2 = biseziona(pol, polB2, dom2, val, toll);
    % Ricorro sulla seconda metà dell'intervallo
    
    D_res = cat(2, D1, D2);
    % Unifico i risultati in una matrice 2xN
    
end


end
