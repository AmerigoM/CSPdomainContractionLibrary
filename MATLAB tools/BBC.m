function dom = BBC(pol, polSimb, dom)
% Implementazione dell'algoritmo BBC
% pol: matrice dei coefficienti del polinomio bi-dimensionale in forma di
%      potenze
% polSimb: polinomio simbolico bi-dimensionale (usare x e y per le variabili)
% dom: vettore 2x2 dei domini delle variabili (espressi per colonne)

% Osservazione: zero non deve appartenere al dominio della seconda
% variabile

x = sym('x','real');
y = sym('y','real');
% definizione variabili simboliche

polB = toBernsteinDuo(pol,dom);
% Trasformazione in Forma di Bernstein

dim = size(polB);
% Dimensione della matrice

% PRIMO ESTREMO PER LA PRIMA VARIABILE

temp = polB(1,:);
g_a1 = [min(temp), max(temp)];
% Calcolo di g(a) per la prima variabile

if(g_a1(1) > 0 && g_a1(2) > 0) % posso contrarre il dominio
    dx = diff(polSimb,x);
    % Calcolo la derivata parziale
    dx = subs(dx,'x',dom(1,1));
    % Sostituisco x con a
    if(strcmp(class(dx),'sym'))
    % Il risultato della sostituzione è simbolico
        polR = sym2poly(dx);
    else
    % Il risultato della sostituzione non è simbolico
        polR = dx;
    end
    for k = 1 : length(polR)
        polRes(k) = polR(length(polR)-(k-1));
    end
    % Converto in forma di coefficienti
    polBB = toBernsteinDuo(polRes,dom);
    % Porto in forma di Bernstein
    g_der1 = [min(polBB), max(polBB)];
    % Calcolo l'intervallo di g'
    
    sol_parziale(1) = min([ g_a1(1)/g_der1(1) g_a1(2)/g_der1(1) g_a1(1)/g_der1(2) g_a1(2)/g_der1(2) ]);
    sol_parziale(2) = max([ g_a1(1)/g_der1(1) g_a1(2)/g_der1(1) g_a1(1)/g_der1(2) g_a1(2)/g_der1(2) ]);
    
    sol_finale = [dom(1,1)-sol_parziale(2), dom(2,1)-sol_parziale(1)];
    % Applicazione della formula di Bernstein-Newton
    
    if(sol_finale(2) >= dom(1,1) && sol_finale(1) <= dom(1,1)) % c'è intersezione
        nuovo_dom(1) = dom(1,1);
        nuovo_dom(2) = sol_finale(2);
    elseif(dom(2,1) >= sol_finale(1) && dom(1,1) <= sol_finale(1)) % c'è intersezione
        nuovo_dom(1) = sol_finale(1);
        nuovo_dom(2) = dom(2,1);
    else
        nuovo_dom(1) = dom(1,1);
        nuovo_dom(2) = dom(2,1);
    end
    
    dom(1,1) = nuovo_dom(1);
    dom(2,1) = nuovo_dom(2);
    % Sostituisco il dominio con la sua contrazione
end

% SECONDO ESTREMO PER LA PRIMA VARIABILE

temp = polB(dim(1),:);
g_b1 = [min(temp), max(temp)];
% Calcolo di g(b) per la prima variabile

if(g_b1(1) > 0 && g_b1(2) > 0) % posso contrarre il dominio
    dx = diff(polSimb,x);
    % Calcolo la derivata parziale
    dx = subs(dx,'x',dom(2,1));
    % Sostituisco x con b
    if(strcmp(class(dx),'sym'))
    % Il risultato è simbolico
        polR = sym2poly(dx);
    else
    % Se il risultato non era simbolico
        polR = dx;
    end
    for k = 1 : length(polR)
        polRes(k) = polR(length(polR)-(k-1));
    end
    % Converto in forma di coefficienti
    polBB = toBernsteinDuo(polRes,dom);
    % Porto in forma di Bernstein
    g_der1 = [min(polBB), max(polBB)];
    % Calcolo l'intervallo di g'
    
    sol_parziale(1) = min([ g_b1(1)/g_der1(1) g_b1(2)/g_der1(1) g_b1(1)/g_der1(2) g_b1(2)/g_der1(2) ]);
    sol_parziale(2) = max([ g_b1(1)/g_der1(1) g_b1(2)/g_der1(1) g_b1(1)/g_der1(2) g_b1(2)/g_der1(2) ]);
    
    sol_finale = [dom(1,1)-sol_parziale(2), dom(2,1)-sol_parziale(1)];
    % Applicazione della formula di Bernstein-Newton
    
    if(sol_finale(2) >= dom(1,1) && sol_finale(1) <= dom(1,1)) % c'è intersezione
        nuovo_dom(1) = dom(1,1);
        nuovo_dom(2) = sol_finale(2);
    elseif(dom(2,1) >= sol_finale(1) && dom(1,1) <= sol_finale(1)) % c'è intersezione
        nuovo_dom(1) = sol_finale(1);
        nuovo_dom(2) = dom(2,1);
    else
        nuovo_dom(1) = dom(1,1);
        nuovo_dom(2) = dom(2,1);
    end
    
    dom(1,1) = nuovo_dom(1);
    dom(2,1) = nuovo_dom(2);
    % Sostituisco il dominio con la sua contrazione
end

% PRIMO ESTREMO PER LA SECONDA VARIABILE

temp = polB(:,1);
g_a2 = [min(temp), max(temp)];
% Calcolo di g(a) per la seconda variabile

if(g_a2(1) > 0 && g_a2(2) > 0) % posso contrarre il dominio
    dy = diff(polSimb,y);
    % Calcolo la derivata parziale
    dy = subs(dy,'y',dom(1,2));
    % Sostituisco y con a
    if(strcmp(class(dy),'sym'))
    % Il risultato è simbolico
        polR = sym2poly(dy);
    else
    % Il risultato non è simbolico
        polR = dy;
    end
    for k = 1 : length(polR)
        polRes(k) = polR(length(polR)-(k-1));
    end
    % Converto in forma di coefficienti
    dom_temp(1,1) = dom(1,2);
    dom_temp(2,1) = dom(2,2);
    dom_temp(1,2) = dom(1,1);
    dom_temp(2,2) = dom(2,1);
    polBB = toBernsteinDuo(polRes,dom_temp);
    % Porto in forma di Bernstein
    g_der1 = [min(polBB), max(polBB)];
    % Calcolo l'intervallo di g'
    
    sol_parziale(1) = min([ g_a2(1)/g_der1(1) g_a2(2)/g_der1(1) g_a2(1)/g_der1(2) g_a2(2)/g_der1(2) ]);
    sol_parziale(2) = max([ g_a2(1)/g_der1(1) g_a2(2)/g_der1(1) g_a2(1)/g_der1(2) g_a2(2)/g_der1(2) ]);
    
    sol_finale = [dom(1,2)-sol_parziale(2), dom(2,2)-sol_parziale(1)];
    % Applicazione della formula di Bernstein-Newton
    
    if(sol_finale(2) >= dom(1,2) && sol_finale(1) <= dom(1,2)) % c'è intersezione
        nuovo_dom(1) = dom(1,2);
        nuovo_dom(2) = sol_finale(2);
    elseif(dom(2,2) >= sol_finale(1) && dom(1,2) <= sol_finale(1)) % c'è intersezione
        nuovo_dom(1) = sol_finale(1);
        nuovo_dom(2) = dom(2,2);
    else
        nuovo_dom(1) = dom(1,2);
        nuovo_dom(2) = dom(2,2);
    end
    
    dom(1,2) = nuovo_dom(1);
    dom(2,2) = nuovo_dom(2);
    % Sostituisco il dominio con la sua contrazione
end

% SECONDO ESTREMO PER LA SECONDA VARIABILE

temp = polB(:,dim(2));
g_b2 = [min(temp), max(temp)];
% Calcolo di g(b) per la seconda variabile

if(g_b2(1) > 0 && g_b2(2) > 0) % posso contrarre il dominio
    dy = diff(polSimb,y);
    % Calcolo la derivata parziale
    dy = subs(dy,'y',dom(2,2));
    % Sostituisco y con b
    if(strcmp(class(dy),'sym'))
    % Se il risultato è simbolico (ossia se sono rimaste entrambe le
    % variabili nella derivazione parziale)
        polR = sym2poly(dy);
    else
    % In questo caso è rimasto un solo valore costante
        polR = dy;
    end
    for k = 1 : length(polR)
        polRes(k) = polR(length(polR)-(k-1));
    end
    % Converto in forma di coefficienti
    dom_temp(1,1) = dom(1,2);
    dom_temp(2,1) = dom(2,2);
    dom_temp(1,2) = dom(1,1);
    dom_temp(2,2) = dom(2,1);
    polBB = toBernsteinDuo(polRes,dom_temp);
    % Porto in forma di Bernstein
    g_der1 = [min(polBB), max(polBB)];
    % Calcolo l'intervallo di g'
    
    sol_parziale(1) = min([ g_b2(1)/g_der1(1) g_b2(2)/g_der1(1) g_b2(1)/g_der1(2) g_b2(2)/g_der1(2) ]);
    sol_parziale(2) = max([ g_b2(1)/g_der1(1) g_b2(2)/g_der1(1) g_b2(1)/g_der1(2) g_b2(2)/g_der1(2) ]);
    
    sol_finale = [dom(1,2)-sol_parziale(2), dom(2,2)-sol_parziale(1)];
    % Applicazione della formula di Bernstein-Newton
    
    if(sol_finale(2) >= dom(1,2) && sol_finale(1) <= dom(1,2)) % c'è intersezione
        nuovo_dom(1) = dom(1,2);
        nuovo_dom(2) = sol_finale(2);
    elseif(dom(2,2) >= sol_finale(1) && dom(1,2) <= sol_finale(1)) % c'è intersezione
        nuovo_dom(1) = sol_finale(1);
        nuovo_dom(2) = dom(2,2);
    else
        nuovo_dom(1) = dom(1,2);
        nuovo_dom(2) = dom(2,2);
    end
    
    dom(1,2) = nuovo_dom(1);
    dom(2,2) = nuovo_dom(2);
    % Sostituisco il dominio con la sua contrazione
end

end

