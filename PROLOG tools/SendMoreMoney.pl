:- use_module(library('clp/bounds')).

puzzle(S,E,N,D,M,O,R,Y):-
    Variabili = [S,E,N,D,M,O,R,Y],
    Variabili in 0..9,
    all_different(Variabili),
    S #>= 1,
    M #>= 1,
    S*1000 + E*100 + N*10 + D + M*1000 + O*100 + R*10 + E #=
    M*10000 + O*1000 + N*100 + E*10 + Y,
    label(Variabili).
