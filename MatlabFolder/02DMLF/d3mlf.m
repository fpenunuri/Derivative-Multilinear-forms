%fr = d3mlf(f,x,y,z,q) computes the derivative multilinear form of rank 3
%for the function f along vectors x, y, z, evaluated at point q.
%f is a dual function D^m --> D^n (or D^(pm) --> D^(pn))
%x, y, z and q are a real arrays of the same dimensions
%x = [x11, ... ,x1m], alternatively
%x = [x11, ... ,x1m; x21, ... ,x2m; ... ;xp1, ... ,xpm]

function fr = d3mlf(f,x,y,z,q)
    fr = 0.5.*(d32mlf(f,x+y,z,q) - d32mlf(f,x,z,q) - d32mlf(f,y,z,q));
end

function fr = d32mlf(f,x,z,q)    
    fr = (d31mlf(f,x+z,q) + d31mlf(f,z-x,q) - 2.*d31mlf(f,z,q))./6;
end

function fr = d31mlf(f,x,q)
    eps1 =  dual3(0,1,0,0);
    fqd = f(q + eps1.*x);
    fr = fqd.f3;
end
