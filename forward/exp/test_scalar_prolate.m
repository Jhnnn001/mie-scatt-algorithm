function test_scalar_prolate
% Independent limits, boundary residual, and analytic derivative checks.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..','spheroid-analytic-forward'));
p = [.3,0,1.5; .7,0,1.1; 1.3,0,.2];
s = scalar_prolate(1,1,1,1,.7,30);
[u,du] = scalar_prolate_eval(s,p);
assert(norm(u-exp(2i*pi*p(:,3))) < 1e-9);
assert(norm(du(:,3)-2i*pi*u) < 1e-8);
% Independent scalar spherical partial waves; the shape differs by 1e-6.
s = scalar_prolate(1,1.12,1,1+1e-6,1,30);
[u,du] = scalar_prolate_eval(s,p);
ref = scalar_sphere(1,1.12,1,1,p);
gap = norm(u-ref)/norm(ref);
assert(gap < 3e-5,'Near-sphere independent reference gap %.3g.',gap);
h = 1e-5;
for dim = 1:3
    pp = p; pm = p; pp(:,dim) = pp(:,dim)+h; pm(:,dim) = pm(:,dim)-h;
    fd = (scalar_prolate_eval(s,pp)-scalar_prolate_eval(s,pm))/(2*h);
    assert(norm(fd-du(:,dim))/norm(du(:)) < 1e-7);
end
fprintf('Scalar reference PASS: independent sphere gap %.3g; boundary %.3g.\n',gap,s.residual);
end

function u = scalar_sphere(lambda,np,nm,a,p)
k = 2*pi*nm/lambda; kp = 2*pi*np/lambda;
ell = (0:30)'; x = k*a; y = kp*a;
jx = sj(ell,x); jy = sj(ell,y); hx = sh(ell,x);
dj = ell/x.*jx-sj(ell+1,x); dh = ell/x.*hx-sh(ell+1,x);
di = ell/y.*jy-sj(ell+1,y);
b = -(k*dj.*jy-kp*jx.*di)./(k*dh.*jy-kp*hx.*di);
r = vecnorm(p,2,2); eta = p(:,3)./r; P = zeros(31,numel(r));
P(1,:) = 1; P(2,:) = eta.';
for l = 1:29
    P(l+2,:) = ((2*l+1)*eta.'.*P(l+1,:)-l*P(l,:))/(l+1);
end
u = exp(1i*k*p(:,3)) + (((2*ell+1).*1i.^ell.*b).'*(sh(ell,k*r.').*P)).';
end

function j = sj(l,z)
j = sqrt(pi./(2*z)).*besselj(l+.5,z);
end

function h = sh(l,z)
h = sqrt(pi./(2*z)).*besselh(l+.5+zeros(size(z)),1,z+zeros(size(l)));
end
