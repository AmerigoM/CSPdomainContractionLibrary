function final = scansionaN(tensore, n_ric, final, A)
    for j = 1 : length(tensore)
        A(n_ric) = j-1;
        if(isnan(tensore(j).value)) % C'è un altra dimensione
            final = scansionaN(tensore(j).next, n_ric + 1, final, A);
        else % Le dimensioni sono finite
            for k = 1 : n_ric % Aggiorno "final" con quello che ho trovato
                if tensore(j).value ~= 0 && A(k) ~= 0 && j-1 > final(k)
                    final(k) = j-1;
                end
            end
        end
    end


end

