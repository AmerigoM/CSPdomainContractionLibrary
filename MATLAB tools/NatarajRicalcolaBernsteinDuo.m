function [B1, B2] = NatarajRicalcolaBernsteinDuo(M, degree, taglio)
% NatarajRicalcolaBernsteinDuo(M, degree, taglio)
% Function che calcola i nuovi polinomi di Bernstein sui nuovi intervalli
% M: matrice dei coefficienti di Bernstein
% degree: grado massimo della variabile di taglio
% taglio: direzione in cui avviene il taglio
%        - 1 : il taglio avviene in corrispondenza del primo intervallo
%        - 2 : il taglio avviene in corrispondenza del secondo intervallo
% B1: matrice dei nuovi coefficienti di Bernstein sul primo tratto
% B2: matrice dei nuovi coefficienti di Bernstein sul secondo tratto

% Controllo sui dati di input
if(taglio ~= 1 && taglio ~= 2)
    fprintf('Errore nei dati di input.\n')
    B1 = NaN;
    B2 = NaN;
    return
end

% Prendo le dimensioni della matrice
dim = size(M);

% Inserisco B^0 nella multimatrice di soluzioni parziali
B_parziali(:,:,1) = M;

% A seconda di come è stata fatta la suddivisione...
if(taglio == 1) % Ho tagliato secondo la prima variabile
    for i = 2 : degree+1 % Indice dei B^i (livello)
        
        for j = 1 : dim(1) % Indice di riga
            for k = 1 : dim(2) % Indice di colonna
                
                if(j < i) % Caso 1 della formula di Nataraj in due demensioni
                    B_parziali(j,k,i) = B_parziali(j,k,i-1);
                else % Caso 2 della formula di Nataraj in due dimensioni
                    B_parziali(j,k,i) = 0.5*B_parziali(j-1,k,i-1) + 0.5*B_parziali(j,k,i-1);
                end
                
            end % End for di colonna
        end % End for di riga
        
    end
    
    % Inserisco la prima soluzione finale in una matrice semplice
    B1 = B_parziali(:,:,degree+1);
    
    % Calcolo B2 per proprietà riflessiva
    for j = 1 : dim(1)
        for k = 1 : dim(2)
            B2temp(dim(1)-(j-1),k) = B_parziali(degree+1,k,(degree+1)-(j-1));
        end
    end
    
    for j = 1 : dim(1)
        for k = 1 : dim(2)
            B2(j,k) = B2temp(dim(1)-(j-1),k);
        end
    end
    
else % Ho tagliato secondo la seconda variabile
    for i = 2 : degree+1 % Indice dei B^i
        
        for j = 1 : dim(1) % Indice di riga
            for k = 1 : dim(2) % Indice di colonna
                
                if(k < i) % Caso 1 della formula di Nataraj in due demensioni
                    B_parziali(j,k,i) = B_parziali(j,k,i-1);
                else % Caso 2 della formula di Nataraj in due dimensioni
                    B_parziali(j,k,i) = 0.5*B_parziali(j,k-1,i-1) + 0.5*B_parziali(j,k,i-1);
                end
                
            end % End for di colonna
        end % End for di riga
        
    end
    
    % Inserisco la prima soluzione finale in una matrice semplice
    B1 = B_parziali(:,:,degree+1);
    
    % Calcolo B2 per proprietà riflessiva
    for j = 1 : dim(2)
        for k = 1 : dim(1)
            B2temp(k,(degree+1)-(j-1)) = B_parziali(k,degree+1,(degree+1)-(j-1));
        end
    end
    
     for j = 1 : dim(1)
        for k = 1 : dim(2)
            B2(j,k) = B2temp(j,dim(2)-(k-1));
        end
    end
    
end

end

