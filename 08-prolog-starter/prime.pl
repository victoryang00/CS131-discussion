% The following functions are potentially useful.

add(A, B, C) :- C is A + B.
lt(A, B) :- A < B.
div(A, B, C) :- C is A / B.
divisible(X, Y) :- div(X, Y, Z), integer(Z).

% Your task is to define a function prime(X) that returns true if X is prime and false if X is not prime.
