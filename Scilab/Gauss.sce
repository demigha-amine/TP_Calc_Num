function [x] = Gauss(T)
   [n,p] = size(T);
   x = zeros(n,p)
   d = T(1,1)
   for i = 2:n
            j = i
            T(i,i) = T(i,i) - ((T(i,j-1) * T(i-1,j)) / d)
            d = T(i,i)
            T(i,j-1) = 0
            printf("%d\t",d)
       
   end
   x = T
endfunction
