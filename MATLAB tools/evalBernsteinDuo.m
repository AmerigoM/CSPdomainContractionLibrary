function res = evalBernsteinDuo(i,j,n1,n2,x,y,a,b,c,d)
% Function per valutare polinomi di Bernstein in due dimensioni in un punto
% n1, n2: gradi massimi delle due variabili
% a,b: dominio della prima variabile
% c,d: dominio della seconda variabile
% x,y: punto di valutazione
% i,j: indici del Bernstein

res1 = nchoosek(n1,i)*( (((b-x)^(n1-i))*((x-a)^i))/((b-a)^n1) );
res2 = nchoosek(n2,j)*( (((d-y)^(n2-j))*((y-c)^j))/((d-c)^n2) );
res = res1*res2;

end