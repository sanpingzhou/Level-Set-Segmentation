function g = Neumann_Bound_Cond(f);

[m, n] = size(f);
g = f;
g([1 m], [1 n]) = g([3 m-2], [3 n-2]);
g([1 m], 2 : end -1) = g([3 m-2], 2 : end - 1);
g(2 : end - 1, [1 n]) = g(2 : end -1, [3 n-2]);

end