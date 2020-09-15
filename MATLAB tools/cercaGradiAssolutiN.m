function final = cercaGradiAssolutiN(tensore, dim)
% C: struttura che implementa una multimatrice
% dim: multi-dimensione del tensore

% Variabile che descrive il numero di volte che ricorro e che segna le
% posizioni all'interno dell'array risultante "final"
n_ric = 1;

% Inizializzaizone di final
final(dim) = 0;

% Inizio scansione
for i = 1 : length(tensore)
    if(~isnan(tensore(i).value)) % Caso monodimensionale
        for k = 1 : length(tensore)
            if(tensore(k).value ~= 0)
                final(i) = k-1;
            end
        end
        return
    else % Caso multidimensionale
        A(n_ric) = i-1; % Variabile utile per i casi di bordo
        final = scansionaN(tensore(i).next, n_ric + 1, final, A);
    end
end
        
        

end

