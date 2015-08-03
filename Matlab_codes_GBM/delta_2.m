function val = delta_2(t1,x,y,h,u_1,u_2)

w_2 = h*B_2(t1,x,y,h,u_1,u_2)/C_2(t1,x,y,h);

if(sign(w_2) >= 0)
    if(abs(w_2)<10^-12)
        val = 0.5;
    elseif(abs(C_2(t1,x,y,h))<10^-12)
        val = 0;
    else
        val = (1/w_2)-(1/(exp(w_2)-1));
    end
else
    if(abs(w_2)<10^-12)
        val = 0.5;
    elseif(abs(C_2(t1,x,y,h))<10^-12)
        val = 1;
    else
        val = (1/w_2)-(1/(exp(w_2)-1));
    end
end