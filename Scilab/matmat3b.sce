function [T]=matmat3b(A,B)
tic()
[m,p] = size(A)
[p,n] = size(B)
T=zeros(m,n)

    for i = 1 : m
        for j = 1 : n
            for k = 1 : p
                T(i, j) = A(i, k )*B(k , j) + T(i, j)
            end
        end
    end
temps=toc()
printf("Time (tic-toc) %f \n",temps);
endfunction





