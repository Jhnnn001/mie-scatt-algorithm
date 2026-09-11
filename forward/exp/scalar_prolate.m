function s = scalar_prolate(lambda,np,nm,a,b,L)
% Experimental scalar Helmholtz reference: U and normal derivative continuous.
% ponytail: axial incidence and real prolate only; not a general public solver.
assert(isreal(np) && np>0 && nm>0 && a>b && b>0 && lambda>0);
s.focal = sqrt(a*a-b*b); s.xi0 = a/s.focal;
s.k = 2*pi*nm/lambda; s.kp = 2*pi*np/lambda;
s.a = a; s.b = b; s.L = L;
co = s.k*s.focal; ci = s.kp*s.focal;
s.eo = sph_eigen('prolate',0,co,L);
s.ei = sph_eigen('prolate',0,ci,L);
[s.rod,chk] = sph_radial_ode(s.eo,co,s.xi0,max(2,s.xi0));
assert(chk.ok,'Outgoing radial ODE failed.');
[eta,w] = gauss_legendre(2*L+24);
[M,rhs] = boundary(s,eta);
weights = sqrt([w;w]); M = M.*weights; rhs = rhs.*weights;
scales = vecnorm(M); c = (M./scales)\rhs; c = c./scales.';
s.bo = c(1:L+1); s.bi = c(L+2:end);
[eta,~] = gauss_legendre(2*L+35);
[M,rhs] = boundary(s,eta);
s.residual = norm(M*c-rhs)/norm(rhs);
assert(s.residual < 1e-7,'Scalar boundary residual %.3g at L=%d.',s.residual,L);
end

function [M,rhs] = boundary(s,eta)
g = struct('eta',eta,'beta',1-eta.^2,'xi',s.xi0,'alpha',s.xi0^2-1);
Wo = sph_angular(s.eo,g); Wi = sph_angular(s.ei,g);
Vo = sph_radial_ode_eval(s.rod,g);
[Vi,ok] = sph_radial(s.ei,1,g,s.kp*s.focal);
assert(ok,'Interior regular radial series failed.');
co = s.k*s.focal;
M = [(Vo.V0.*Wo.W0).', -(Vi.V0.*Wi.W0).'; ...
    (Vo.V1.*Wo.W0).'/sqrt(g.alpha)/co, ...
    -(Vi.V1.*Wi.W0).'/sqrt(g.alpha)/co];
inc = exp(1i*co*s.xi0*eta);
rhs = [-inc; -1i*eta.*inc];
end
