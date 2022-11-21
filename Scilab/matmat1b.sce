function [T]=matmat1b(A,B)
tic()
[m,p] = size(A)
[p,n] = size(B)
T=zeros(m,n)

        
        for i = 1 : m            
               T(i, :) = A(i, :)*B + T(i, :);
                
        end 
temps=toc()
printf("Time (tic-toc) %f \n",temps);

endfunction




