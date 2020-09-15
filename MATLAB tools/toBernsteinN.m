function B_ff = toBernsteinN(tensore, dom, ass)
% C: struttura tensoriale rappresentante il polinomio
% dom: matrice 2*N contenente i domini delle singole variabili
% ass: vettore dei gradi assoluti delle singole variabili
% B_f: Struttura contenente la forma implicita dei coefficienti di
%      Bernstein

% A che livello di ricorsione sono arrivato
n_ric = 1;

% Multi-matrice risultante
B_f.f = 'Forma di Bernstein Implicita';

% Variabile globale
global contatore;
contatore = 1;

% Struttura finale
B_f.next = B_f;

% Inizio scansione
for i = 1 : length(tensore)
    if(~isnan(tensore(i).value)) % Caso monodimensionale
        for k = 1 : length(tensore)
            vett1dim(k) = tensore(k).value;
        end
        B_ff = toBernsteinMono(vett1dim);
        return
    else % Caso multidimensionale
        rel(n_ric) = i-1; % Vettore dei gradi relativi
        B_f.next = monomioN(tensore(i).next, dom, ass, rel, n_ric + 1, B_f.next);
        
    end
end

% Pulisco il risultato
B_f = B_f.next;

counter = 1;
for k = 1 : contatore - 1
    if(isfield(B_f,(['coeff' num2str(k)])))
        B_ff.(['coeff' num2str(counter)]) = B_f.(['coeff' num2str(k)]);
    end
    
    if(~isnan( B_f.(['f' num2str(k)]) ))
        B_ff.(['f' num2str(counter)]) = B_f.(['f' num2str(k)]);
        counter = counter + 1;
    end
end

end

