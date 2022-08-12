%fr = d2mlf(f,x,y,q) computes the derivative multilinear form of rank 2
%for the function f along vectors x and y, evaluated at point q.
%f is a dual function D^m --> D^n (or D^(pm) --> D^(pn))
%x, y and q are a real arrays of the same dimensions
%x = [x11, ... ,x1m], alternatively
%x = [x11, ... ,x1m; x21, ... ,x2m; ... ;xp1, ... ,xpm]

function fr = d2mlf(f,x,y,q)
    fr = 0.5.*(d21mlf(f,x+y,q) - d21mlf(f,x,q) - d21mlf(f,y,q));
end

function fr = d21mlf(f,v,q)
    eps1 =  dual2(0,1,0);
    fqd = f(q + eps1.*v);
    fr = fqd.f2;
end
