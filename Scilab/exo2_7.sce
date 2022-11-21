A=rand(20,20)
U= triu(A)
L= tril(A)


xex=rand(20,1)

b=U * xex
bb=L * xex

x_res=usolve(U,b)

err_av_u=norm(xex - x_res)/norm(xex)
err_ar_l=norm(b-A*x_res)/(norm(A)*norm(x_res))


pr_u=cond(U)*err_ar_u


x_r=usolve(U,bb)

err_av_u=norm(xex - x_r)/norm(xex)
err_ar_l=norm(bb-A*x_r)/(norm(A)*norm(x_r))


pr_l=cond(U)*err_ar_l
