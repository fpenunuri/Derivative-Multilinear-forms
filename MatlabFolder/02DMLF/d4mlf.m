%fr = d4mlf(f,x,y,z,w,q) computes the derivative multilinear form of rank 4
%for the function f along vectors x, y, z, w, evaluated at point q.
%f is a dual function D^m --> D^n (or D^(pm) --> D^(pn))
%x, y, z, w and q are a real arrays of the same dimensions
%x = [x11, ... ,x1m], alternatively
%x = [x11, ... ,x1m; x21, ... ,x2m; ... ;xp1, ... ,xpm]

function fr = d4mlf(f,x,y,z,w,q)
    fr = 0.5.*(d43mlf(f,x+y,z,w,q) - d43mlf(f,x,z,w,q) - ...
        d43mlf(f,y,z,w,q));
end

function fr = d43mlf(f,x,z,w,q)
    fr = 0.5.*(d42mlf(f,x,z+w,q) - d42mlf(f,x,z,q) - d42mlf(f,x,w,q));
end

function fr = d42mlf(f,x,z,q)    
    d4xzp = d41mlf(f, x + z, q);
    d4xzm = d41mlf(f, x - z, q);
    fr = (d4xzp + d4xzm - 2.*d41mlf(f,x,q) - 2.*d41mlf(f,z,q))./12;
end

function fr = d41mlf(f,x,q)
    eps1 =  dual4(0,1,0,0,0);
    fqd = f(q + eps1.*x);
    fr = fqd.f4;
end
