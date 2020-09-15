function B = monomioN(tensore, dom, ass, rel, n_ric, B)
% Function che ricorre su ogni lista della struttura
% tensore: struttura multimatrice
% dom: dominio di definizione delle variabili
% ass: vettore dei gradi assoluti
% rel: indice di ricorsione (permette di tenere in conto il grado relativo)
% n_ric: indice di rel

global contatore;

for j = 1 : length(tensore)
    rel(n_ric) = j-1;
    if(isnan(tensore(j).value)) % C'è un'altra dimensione
        B = monomioN(tensore(j).next, dom, ass, rel, n_ric + 1, B);
    else
        if(~isnan(tensore(j).value))
            if(tensore(j).value ~= 0)
                for k = 1 : n_ric
                    eval(['B' num2str(k) '= monomio(rel(k),ass(k),dom(1,k),dom(2,k));'])
                    % I Bi sono i vettori ottenuti dai monomi da
                    % moltiplicare tensorialmente fra di loro
                end
                
                B.( ['coeff' num2str(contatore)] ) = tensore(j).value;
                for k = 1 : n_ric
                    B.( ['f' num2str(contatore)] ) = eval(['B' num2str(k)]);
                    contatore = contatore + 1;
                end
                
            else
               B.( ['f' num2str(contatore)] ) = NaN;
               contatore = contatore + 1;
            end
            
        else
             B.( ['f' num2str(contatore)] ) = NaN;
             contatore = contatore + 1;
        end
    end
    
end


B.( ['f' num2str(contatore)] ) = NaN;
contatore = contatore + 1;

end
