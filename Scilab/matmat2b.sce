function [T]=matmat2b(A,B)
tic()
[m,p] = size(A)
[p,n] = size(B)
T=zeros(m,n)

  
        
        for i = 1 : m
            for j = 1 : n
                T(i, j) = A(i, :)*B(:, j) + T(i, j);
                
            end 
        end
temps=toc()
printf("Time (tic-toc) %f \n",temps);
endfunction
