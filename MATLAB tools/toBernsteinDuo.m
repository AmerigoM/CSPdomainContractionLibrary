function B_f = toBernsteinDuo(C, dom)
% Function per la conversione di un polinomio in due
% variabili da Forma di Potenze a Forma di Bernstein:
% B_f = toBernsteinDuo(C, dom)
% C: matrice dei coefficienti del polinomio in due variabili
%    (dove il termine C(i,k) è il coefficiente di x^(i-1)*y^(k-1))
% dom: matrice 2*n contenente i domini delle singole variabili. Gli
%      estremi del primo dominio, ad esempio, sono in posizione dom(1,1)
%      e dom(2,1)
% B_f: matrice dei coefficienti del polinomio in base di Bernstein

g_ass = [0 0];
% Vettore dei gradi assoluti delle due variabili

dimM = size(C);
% Calcolo le dimenzioni della matrice in input

% Calcolo i gradi assoluti delle singole variabili
for i = 1 : dimM(1)
    for j = 1 : dimM(2)
        if(g_ass(1) < i && C(i,j) ~= 0)
            g_ass(1) = i-1;
        end
        
        if(g_ass(2) < j && C(i,j) ~= 0)
            g_ass(2) = j-1;
        end
    end
end

B_f = 0;
% Inizializzo la matrice finale dei coefficienti di Bernstein del
% polinomio in due variabili

for i = 1 : dimM(1)
    for j = 1 : dimM(2)
        if(C(i,j) ~= 0 && (i ~= 1 || j ~= 1) )
        % Se l'elemento nella matrice è non nullo e non mi trovo nè sulla
        % prima riga né sulla prima colonna (non sono monomi insomma)
            B1 = monomio(i-1,g_ass(1),dom(1,1),dom(2,1));
            % Calcolo il vettore di coeff di Bernstein nella prima variabile
            B2 = monomio(j-1,g_ass(2),dom(1,2),dom(2,2));
            % Calcolo il vettore di coeff di Bernstein nell'altra variabile
            B = C(i,j)*B1'*B2;
            % Calcolo i coeff di Bernstein del monomio in più variabili
            B_f = B_f + B;
        elseif(C(i,j) ~= 0 && i == 1 && j ~= 1)
        % Caso in cui è presente solo la seconda variabile
            B2 = monomio(j-1,g_ass(2),dom(1,2),dom(2,2));
            B_f = B_f + C(i,j)*ones(g_ass(1),1)*B2;
        elseif(C(i,j) ~= 0 && j == 1 && i ~= 1)
        % Caso in cui è presente solo la prima variabile
        	B1 = monomio(i-1,g_ass(1),dom(1,1),dom(2,1));
            B_f = B_f + C(i,j)*B1*ones(1,g_ass(2));
        elseif(C(i,j) ~= 0 && j == 1 && i == 1)
        % Caso di bordo
            B_f = B_f + C(1,1)*ones(g_ass(1)+1,1)*ones(1,g_ass(2)+1);
        end
    end
end


end

