lrevaux([], A, A).
lrevaux([X | Y], A, R) :- lrevaux(Y, [X | A], R).
lrev(X, R) :- lrevaux(X, [], R).
