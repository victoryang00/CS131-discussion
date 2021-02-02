/* typeof(X, T, Env): Type of X is T in */
/*     environment Env. */

typeof(X, T, Env) :- defn(X, T, Env).
typeof(X, int, _) :- integer(X).

typeof([], [_], _).
typeof([E | R], [T], Env) :- typeof(E, T, Env), typeof(R, [T], Env).
typeof(L + R, int, Env) :- typeof(L, int, Env), typeof(R, int, Env).
typeof(lambda(X, E), T1->T2, Env) :- typeof(E, T2, [def(X, T1) | Env]).
typeof(L // R, [T], Env) :- typeof(L, [T], Env), typeof(R, [T], Env). /* Think of this as performing some type-preserving operations on lists, such as concatenation. */
typeof(L << R, T2, Env) :- typeof(R, T1, Env), typeof(L, T1->T2, Env). /* Think of this as function application. */
typeof(cast(L, R), T1 -> T2, Env) :- typeof(L, T1, Env), typeof(R, T2, Env). /*Think of this as a transformation from type of L to type of R, such as a cast.*/
