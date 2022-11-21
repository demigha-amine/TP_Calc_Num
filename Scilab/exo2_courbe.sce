cl=zeros(19,1)
cu=zeros(19,1)
for i=2:40
    A=rand(i,i)
    U= triu(A)
    L= tril(A)
    
    
    xex=rand(i,1)
    
    b=U * xex
    bb=L * xex
    
    x_res=usolve(U,b)
    
    err_av_u=norm(xex - x_res)/norm(xex)
    err_ar_u=norm(b-A*x_res)/(norm(A)*norm(x_res))
    
    
    pr_u=cond(U)*err_ar_u
    cu(i)= pr_l - err_av_u
    
    x_r=usolve(U,bb)
    
    err_av_l=norm(xex - x_r)/norm(xex)
    err_ar_l=norm(bb-A*x_r)/(norm(A)*norm(x_r))
    
    
    pr_l=cond(U)*err_ar_l
    cl(i)= pr_l - err_av_l
end
