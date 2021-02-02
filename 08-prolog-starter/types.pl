/* typeof(X, T, Env): Type of X is T in */
/*     environment Env. */

typeof(X, T, Env) :- defn(X, T, Env).
typeof(X, int, _) :- integer(X).
