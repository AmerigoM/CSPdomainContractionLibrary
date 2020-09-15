function [polRes, polResB] = ricalcolaBernstein(pol, min, max)
% Funzione che effettua la traformazione affine [min, max] -> [0,1]
% per un polinomio pol e ricalcola i coefficienti di Bernstein
% [polRes, polResB] = ricalcolaBernstein(pol, min, max)
% pol: vettore di coefficienti del polinomio (in forma di potenze) dal
%      meno significativo al più significativo
% min, max: estremi del dominio
% polRes: polinomio in forma di potenze risultante, espresso dal
%         coefficiente meno significativo al più significativo
% polResB: polinomio polRes in base di Bernstein ricalcolato

A = [0 1; 1 1];
b = [min, max];
sol = b/A;
% Risolvo il sistema linare per compiere la trasformazione affine

grado = length(pol) - 1;
% grado del polinomio

symPol = sym(pol(1));
for i = 1 : grado
    symPol = symPol + sym(pol(i+1))*'x'^sym(i);
end
% Creo il polinomio simbolico

symPol = simplify(symPol);
% Polinomio in forma simbolica risultante
toChange = sym(sol(1))*'x' + sym(sol(2));
% Polinomio da sostituire ad x per completare la trasformazione affine
symPol = subs(symPol, 'x', toChange);
symPol = simplify(symPol);
% Faccio la sostituzione e semplifico
polR = sym2poly(symPol);
% Converto in forma di coefficienti e restituisco

for k = 1 : length(polR)
   polRes(k) = polR(length(polR)-(k-1));
end

polResB = toBernsteinMono(polRes);
% Converto il polniomio calcolato in forma di Bernstein


end

