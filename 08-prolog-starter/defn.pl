/* defn(I,D,Z): I is defined as D in Z. */
/* Z is of the form [def(a,b), ...]. */
defn(I, T, [def(I, T) | _]).
defn(I, T, [def(I1, _) | R]) :- dif(I, I1), defn(I, T, R).
typeof(X, T, Env) :- defn(X, T, Env).

/* 

| ?- defn(c, d, [def(a,b), def(c,d)]).

yes

| ?- defn(a, X, [def(a,b), def(c,d)]).

X = b 

| ?- defn(X, d, [def(a,b), def(c,d)]).

X = c 
*/
