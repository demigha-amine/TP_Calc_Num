A=rand(3,3)


xex=rand(3,1)

b=A*xex

x_res=A\b

err_av=norm(xex-x_res)/norm(xex)
err_ar=norm(b-A*x_res)/(norm(A)*norm(x_res))
pr=cond(A)*err_ar

