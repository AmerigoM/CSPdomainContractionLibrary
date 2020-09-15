function B_f = toBernsteinDuoTest(C, dom)
% Function per la conversione di un polinomio in due
% variabili da Forma di Potenze a Forma di Bernstein:
% B_f = toBernsteinMono(C, dom)
% C: matrice dei coefficienti del polinomio in due variabili
%    (dove il termine C(i,k) è il coefficiente di x^(i-1)*y^(k-1))
% dom: matrice 2*n contenente i domini delle singole variabili
% B_f: matrice dei coefficienti del polinomio in base di Bernstein

g_ass = [0 0];
% Vettore dei gradi assoluti delle due variabili

% Calcolo i gradi assoluti delle singole variabili
for i = 1 : length(C)
    for j = 1 : length(C)
        if(g_ass(1) < i && C(i,j) ~= 0 && i ~= 0)
            g_ass(1) = i-1;
        end
        
        if(g_ass(2) < j && C(i,j) ~= 0 && j ~= 0)
            g_ass(2) = j-1;
        end
    end
end

B_f = 0;
% Inizializzo la matrice finale dei coefficienti di Bernstein del
% polinomio in due variabili

for i = 1 : length(C)
    for j = 1 : length(C)
        if(C(i,j) ~= 0)
        % Se l'elemento nella matrice è non nullo e non mi trovo nè sulla
        % prima riga né sulla prima colonna (non sono monomi insomma)
            B1 = monomio(i-1,g_ass(1),dom(1,1),dom(2,1));
            % Calcolo il vettore di coeff di Bernstein nella prima variabile
            B2 = monomio(j-1,g_ass(2),dom(1,2),dom(2,2));
            % Calcolo il vettore di coeff di Bernstein nell'altra variabile
            B = C(i,j)*B1'*B2;
            % Calcolo i coeff di Bernstein del monomio in più variabili
            B_f = B_f + B;
        end
    end
end


end