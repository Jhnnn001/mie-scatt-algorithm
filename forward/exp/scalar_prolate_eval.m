function [u,du,us] = scalar_prolate_eval(s,p)
% Total scalar field, Cartesian gradient, and scattered field outside.
assert(size(p,2)==3 && all(isfinite(p),'all'));
u = complex(zeros(size(p,1),1)); du = complex(zeros(size(p))); us = u;
for start = 1:1000:size(p,1)
    ix = start:min(start+999,size(p,1));
    g = sph_coords('prolate',p(ix,1)/s.focal,p(ix,2)/s.focal,p(ix,3)/s.focal);
    outer = g.xi >= s.xi0;
    for side = [false,true]
        jj = find(outer==side);
        if isempty(jj), continue; end
        q = sph_coords('prolate',p(ix(jj),1)/s.focal,p(ix(jj),2)/s.focal,p(ix(jj),3)/s.focal);
        if side
            W = sph_angular(s.eo,q); V = sph_radial_ode_eval(s.rod,q); c = s.bo;
        else
            W = sph_angular(s.ei,q); [V,ok] = sph_radial(s.ei,1,q,s.kp*s.focal); c = s.bi;
            assert(ok,'Interior evaluation radial series failed.');
        end
        value = (c.'*(W.W0.*V.V0)).';
        % Metric factors cancel square roots in the regular weighted derivatives.
        dxi = (c.'*(W.W0.*V.V1)).'./(s.focal*sqrt(q.D));
        deta = (c.'*(W.W1.*V.V0)).'./(s.focal*sqrt(q.D));
        grad = dxi.*q.e_xi + deta.*q.e_eta;
        ids = ix(jj);
        if side
            us(ids) = value;
            inc = exp(1i*s.k*p(ids,3)); value = value+inc;
            grad(:,3) = grad(:,3)+1i*s.k*inc;
        end
        u(ids) = value; du(ids,:) = grad;
    end
end
end
