function A = scalar_far_amplitude(s,eta)
% Exact outgoing r*exp(-i*k*r) limit, with no large common phase evaluation.
shape = size(eta); eta = eta(:);
r = s.rod;
c = s.focal*exp(-1i*r.c*r.xi0)*exp(-r.qtilde_xi0(:)) .* ...
    r.asymptotic_coefficients(1,:).' ./ ...
    (r.asymptotic_argument_scale*r.asymptotic_U_start(:));
A = complex(zeros(size(eta)));
for first = 1:2000:numel(eta)
    ix = first:min(first+1999,numel(eta));
    W = sph_angular(s.eo,struct('eta',eta(ix),'beta',1-eta(ix).^2));
    A(ix) = ((s.bo.*c).'*W.W0).';
end
A = reshape(A,shape);
end
