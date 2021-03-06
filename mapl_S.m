function S = mapl_S(n,m)

S = 2*(-1)^n*pi^(3.5)*(kronDel(n-m)*3*(2*n^2+2*n+1) + ...
    kronDel(n+2-m)*(6+4*n)*sqrt(factorial(m)/factorial(n)) + ...
    kronDel(n+4-m)*sqrt(factorial(m)/factorial(n)) + ...
    kronDel(m+2-n)*(6+4*m)*sqrt(factorial(n)/factorial(m)) + ...
    kronDel(m+4-n)*sqrt(factorial(n)/factorial(m)));

end
