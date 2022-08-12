%fr = d1mlf(f,v,q) computes the derivative multilinear form of rank 1
%for the function f along vector v evaluated at point q.
%f is a dual function D^m --> D^n (or D^(pm) --> D^(pn))
%v and q are a real arrays of the same dimensions
%v = [v11, ... ,v1m], alternatively
%v = [v11, ... ,v1m; v21, ... ,v2m; ... ;vp1, ... ,vpm]
function fr = d1mlf(f,v,q)
    eps1 =  dual1(0,1);
    fqd = f(q + eps1.*v);
    fr = fqd.f1;
end
