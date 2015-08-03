function val = inter(u,i,j,c)

if(c==1)
    val = (u(i+1,j) + u(i,j))/2;
else
    val = (u(i,j) + u(i,j+1))/2;
end