function [rhs,i] = richardson(A,b,n,tol)
   [m,p] = size(A)
   rhs = ones(1,m)
   I = eye(m,p) 
   i = 0
   lambda_min = min(spec(A))
   lambda_max = max(spec(A))
   alpha = 2 / (lambda_min + lambda_max)
   while (norm(rhs) > tol & (i < n))
       rhs = rhs * (I - alpha * A) + (alpha * b)
        i = i + 1
   end
  
endfunction
