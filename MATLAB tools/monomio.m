function b = monomio(g_rel, g_ass, minimo, massimo)
% Function che calcola i coefficienti di Bernstein di un monomio con una
% variabile, avente un certo grado assoluto e un certo grado relativo.
% b = monomio(g_rel, g_ass, minimo, massimo)
% b: vettore contenente i coefficienti di Bernstein del monomio
% g_rel: grado relativo del monomio
% g_ass: grado assoluto del monomio
% minimo: estremo inferiore del dominio
% massimo: estremo superiore del dominio

if(g_rel == g_ass) % PRIMO CASO
    for k = 0 : g_ass
        b(k+1) = minimo^(g_rel - k) * massimo^(k);
    end
else % SECONDO CASO
    for k = 0 : g_ass
    % Calcolo tutti i b(k)
        r = g_ass - g_rel;
        j = max([0,k-r]);
        jj = min([k,g_rel]);
        b(k+1) = 0;
        % Setto un po' di parametri per poter applicare la formula in una
        % maniera più immediata
        
        for i = j : jj
        % Applico la formula del SECONDO CASO
           b(k+1) = b(k+1) + ((nchoosek(r, k-i)*nchoosek(g_rel, i))/nchoosek(g_rel+r, k))*(minimo^(g_rel-i)*massimo^i);
        end
        
    end
end

end