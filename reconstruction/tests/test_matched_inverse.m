function test_matched_inverse
% Incorrect FFT scale, interpolation transpose or missing outer update fails here.
assert(exist('odt_operator','file')==2,'Missing matched Ewald operator.');
rng(41); sz=[8,10,6]; dx=.1;
q=(rand(91,3)-.5)*18;
op=odt_operator(q,sz,dx,2,linspace(.1,1,91).');
x=randn(sz)+1i*randn(sz); v=randn(91,1)+1i*randn(91,1);
ax=op.forward(x); ah=op.adjoint(v);
assert(abs(ax'*v-x(:)'*ah(:))/max(norm(ax)*norm(v),eps)<1e-12);
% Pruning unused frequencies must equal the unpruned complex 3-D FFT.
ps=2*sz; padded=complex(zeros(ps)); offset=(ps-sz)/2;
padded(offset(1)+(1:sz(1)),offset(2)+(1:sz(2)),offset(3)+(1:sz(3)))=x;
ft=fftshift(fftn(ifftshift(padded)))*dx^3;
f=arrayfun(@(n) (-n/2:n/2-1)*2*pi/(n*dx),ps,'UniformOutput',false);
expected=interpn(f{1},f{2},f{3},ft,q(:,1),q(:,2),q(:,3),'linear');
assert(norm(ax./op.sample_scale-expected)/norm(expected)<1e-12);
point=zeros(sz); point(5,6,4)=.03;
assert(max(abs(op.forward(point)./op.sample_scale-.03*dx^3))<1e-13);
% Accepted endpoint roundoff must use the same clamping in the crop and stencil.
edge=odt_operator([-pi/dx-1e-12,0,0],[8,8,8],dx,1);
edge_point=zeros(8,8,8); edge_point(5,5,5)=.03;
assert(abs(edge.forward(edge_point)/edge.sample_scale-.03*dx^3)<1e-13);
% Internal interpolation padding must not change the direct reconstruction grid.
native=odt_operator(q,sz,dx,1,linspace(.1,1,91).');
physical=randn(91,1)+1i*randn(91,1);
native_direct=native.direct(physical.*native.sample_scale);
padded_direct=op.direct(physical.*op.sample_scale);
assert(norm(native_direct(:)-padded_direct(:))/norm(native_direct(:))<1e-12, ...
    'Direct gridding changed with internal interpolation padding.');
nm=1.335381534; n=nm*sqrt(1+.01*rand(sz));
[b,~,u0,g]=born_rytov_voxel(.532,n,nm,dx,[0,.12],[0,.7],1,1);
d=odt_prepare(u0+b,u0,g,'born');
chi=n.^2/nm^2-1;
assert(norm(d.operator.forward(chi)-d.y)/norm(d.y)<1e-11);
% Independent, fully sampled and hand-solvable two-level TV problem.
sz=[8,8,8]; f=(-4:3)*2*pi/(8*dx); [QX,QY,QZ]=ndgrid(f,f,f);
fullop=odt_operator([QX(:),QY(:),QZ(:)],sz,dx,1);
x=.03*ones(sz); x(5:8,:,:)=.09;
full=struct('operator',fullop,'y',fullop.forward(x),'n_m',nm);
[~,direct]=odt_reconstruct(full,'direct');
assert(norm(direct(:)-x(:))/norm(x(:))<1e-12);
[~,gp,igp]=odt_reconstruct(full,'gp',struct('max_iter',5));
assert(norm(gp(:)-x(:))/norm(x(:))<1e-10 && igp.converged);
[~,tv,itv]=odt_reconstruct(full,'tv',struct('alpha',.008, ...
    'outer_iter',1,'max_iter',1200,'tol',1e-7,'cg_tol',1e-10));
expected=.034*ones(sz); expected(5:8,:,:)=.086;
assert(max(abs(tv(:)-expected(:)))<3e-6 && itv.inner_converged);
% Scaling data objective, TV and both penalties leaves ADMM iterates unchanged.
scaled=full; factor=.25;
scaled.operator.forward=@(v) sqrt(factor)*fullop.forward(v);
scaled.operator.adjoint=@(v) sqrt(factor)*fullop.adjoint(v);
scaled.operator.direct=@(v) fullop.direct(v/sqrt(factor));
scaled.operator.normal_diagonal=factor*fullop.normal_diagonal;
scaled.y=sqrt(factor)*full.y;
[~,first]=odt_reconstruct(full,'tv',struct('alpha',.008,'outer_iter',1,'max_iter',20,'tol',1e-10));
[~,second]=odt_reconstruct(scaled,'tv',struct('alpha',factor*.008, ...
    'rho',factor*.1,'rho_positive',factor*.1,'outer_iter',1,'max_iter',20,'tol',1e-10));
assert(norm(first(:)-second(:))/norm(first(:))<1e-10);
[~,reinject,iri]=odt_reconstruct(full,'tv',struct('alpha',.008, ...
    'outer_iter',3,'max_iter',1200,'tol',1e-7,'cg_tol',1e-10));
assert(norm(reinject(:)-x(:))<norm(tv(:)-x(:))*.01);
assert(iri.outer_iterations>=2 && iri.reinjection_norm>0);
% Nonuniform GP must agree with an independent real pseudoinverse projection.
tiny=odt_operator((rand(22,3)-.5)*40,[4,4,4],dx,2);
matrix=zeros(44,64);
for j=1:64
    e=zeros(4,4,4); e(j)=1; column=tiny.forward(e);
    matrix(:,j)=[real(column);imag(column)];
end
truth=rand(4,4,4); truth(truth<.6)=0;
td=struct('operator',tiny,'y',tiny.forward(truth),'n_m',nm);
start=max(tiny.direct(td.y),0); residual=td.y-tiny.forward(start);
expected=max(start(:)+pinv(matrix)*[real(residual);imag(residual)],0);
[~,one]=odt_reconstruct(td,'gp',struct('max_iter',1,'cg_max_iter',300,'cg_tol',1e-11));
assert(norm(one(:)-expected)/norm(expected)<1e-6);
[~,~,short]=odt_reconstruct(full,'tv',struct('max_iter',1,'outer_iter',1,'tol',1e-12));
assert(~short.inner_converged);
% A large splitting penalty must not falsely certify the unchanged input.
[~,~,stalled]=odt_reconstruct(full,'tv',struct('alpha',.008,'rho',1000, ...
    'rho_positive',1000,'outer_iter',1,'max_iter',10));
assert(~stalled.inner_converged,'ADMM residuals alone hid a nonstationary solution.');
full.y(:)=0;
for method={'direct','gp','tv'}
    [nz,xz,iz]=odt_reconstruct(full,method{1});
    assert(all(xz(:)==0) && all(nz(:)==nm) && iz.converged);
end
fprintf('Matched adjoint/forward, analytic TV and Lim reinjection: PASS\n');
end
