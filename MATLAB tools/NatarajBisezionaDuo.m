function D_res = NatarajBisezionaDuo(polB, grado1, grado2, dom1, dom2, val, toll)
% Function che implementa la bisezione del dominio
% polB : polinomio in Forma di Bernstein
% grado1: grado massimo della prima variabile
% grado2: grado massimo della seconda variabile
% dom1: dominio della prima variabile
% dom2: dominio della seconda variabile
% val: valore del vincolo
% toll: tolleranza assegnata

massimo = max(max(polB));
minimo = min(min(polB));

if(massimo > val && minimo > val)
    D_res(1,1) = NaN;
    D_res(2,1) = NaN;
    D_res(3,1) = NaN;
    D_res(4,1) = NaN;
elseif(massimo < val+toll && minimo < val+toll)
    D_res = dom1;
    D_res = [D_res; dom2];
else
    if( abs(dom1(2,1) - dom1(1,1)) >= abs(dom2(2,1) - dom2(1,1)) )
    % Taglio secondo la prima variabile
        middle = ((dom1(1,1) + dom1(2,1))/2);
        domA(1,1) = dom1(1,1);
        domA(2,1) = middle;
        domB(1,1) = middle;
        domB(2,1) = dom1(2,1);
        
        [polB1, polB2] = NatarajRicalcolaBernsteinDuo(polB,grado1,1);
        D1 = NatarajBisezionaDuo(polB1, grado1, grado2, domA, dom2, val, toll);
        D2 = NatarajBisezionaDuo(polB2, grado1, grado2, domB, dom2, val, toll);
        
        % Concateno i risultati ottenuti
        D_res = cat(2, D1, D2);
    else
        % Taglio secondo la seconda variabile
        middle = ((dom2(1,1) + dom2(2,1))/2);
        domA(1,1) = dom2(1,1);
        domA(2,1) = middle;
        domB(1,1) = middle;
        domB(2,1) = dom2(2,1);
        
        [polB1, polB2] = NatarajRicalcolaBernsteinDuo(polB,grado2,2);
        D1 = NatarajBisezionaDuo(polB1, grado1, grado2, dom1, domA, val, toll);
        D2 = NatarajBisezionaDuo(polB2, grado1, grado2, dom1, domB, val, toll);
        
        % Concateno i risultati ottenuti
        D_res = cat(2, D1, D2);
    end
    
end

end

