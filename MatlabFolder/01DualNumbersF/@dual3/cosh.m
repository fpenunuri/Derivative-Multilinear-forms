%COSH   Hyperbolic cosine.
% The elements of X can be dual3 numbers
function fr = cosh(g) 
  g0 = g.f0;
  g1 = g.f1;
  g2 = g.f2;
  g3 = g.f3;
  
  fr0 = cosh(g0);
  fr1 = sinh(g0).*g1;
  fr2 = g2.*sinh(g0) + cosh(g0).*(g1.^2);
  fr3 = 3.*g1.*g2.*cosh(g0) + g1.^3 .* sinh(g0) + g3.*sinh(g0);
  
  fr = dual3(fr0,fr1,fr2,fr3);
end
