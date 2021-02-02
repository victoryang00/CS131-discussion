% Students can assume the following are defined.
% Potentially useful functions:

add(A, B, C) :- C is A + B.
lt(A, B) :- A < B.
div(A, B, C) :- C is A / B.
divisible(X, Y) :- div(X, Y, Z), integer(Z).

% We would like to write a predicate, composite(X), that checks if a number
% X is composite (that is, whether it can be represented as Y * Z for Y > 1,
% Z > 1).
%
% To do this, we create a helper predicate, composite(X, Y), that checks if
% X can be divided by Y, or any number greater than Y but smaller than X / 2.
composite(X, Y) :- Y > 1, divisible(X, Y).
composite(X, Y) :- lt(Y, X / 2), composite(X, Y+1).

% Now, our task is much simpler. Define composite(X) using the two argument
% version of composite.
composite(X) :- X > 2, composite(X, 2).
prime(X) :- not(composite(X)).
