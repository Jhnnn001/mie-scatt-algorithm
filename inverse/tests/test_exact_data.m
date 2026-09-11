function test_exact_data
% A wrong far-to-plane scale or polarization projection breaks this comparison.
assert(exist('odt_exact_data','file')==2,'Missing independent exact detector data.');
assert(exist('odt_predict','file')==2,'Missing physical detector prediction.');
m=odt_exact_data(64,.2,1);
assert(isequal(size(m.U),[64,64]) && all(isfinite(m.U),'all'));
assert(abs(m.parameters.n_p-1.33567771866)<1e-12);
assert(norm(m.U(:)-m.u0(:))>0 && m.geometry.NA_det==1);
d=odt_prepare(m.U,m.u0,m.geometry,'rytov');
assert(all(isfinite(d.y)) && max(d.phase_order_error)<1e-6);
% Same detector grid and a genuinely tilted, nonzero-azimuth archived case.
tilted=odt_exact_data(512,.1,21);
dk=2*pi/(512*.1);
co_norm2=sum(abs(tilted.co_spectrum).^2)*dk^2;
scalar_norm2=sum(abs(tilted.scalar_spectrum).^2)*dk^2;
assert(abs(co_norm2/tilted.angles.co_norm2-1)<1e-11);
assert(abs(scalar_norm2/tilted.angles.scalar_norm2-1)<1e-11);
% An independent weak-contrast Born amplitude also checks the complex phase.
gg=tilted.geometry; [kx,ky]=ndgrid(gg.kx,gg.ky); pp=gg.pupil;
q=[kx(pp),ky(pp),gg.kz(pp)]-gg.k_incident;
Q=sqrt(25*(q(:,1).^2+q(:,2).^2)+9*q(:,3).^2);
h=ones(size(Q)); nz=abs(Q)>1e-5;
h(nz)=3*(sin(Q(nz))-Q(nz).*cos(Q(nz)))./Q(nz).^3;
p=tilted.parameters;
born=1i*gg.k0^2*(p.n_p^2-p.n_m^2)*100*pi*h.*exp(1i*gg.kz(pp)*5)./(2*gg.kz(pp));
assert(norm(born-tilted.co_spectrum)/norm(tilted.co_spectrum)<.02);
nm=1.335381534; n=nm*ones(16,16,16); n(7:10,7:10,7:10)=nm+.0003;
[b,r,~,g]=born_rytov_voxel(.532,n,nm,.1,[0,.2],[0,.7],1,2);
[pb,pr]=odt_predict(n.^2/nm^2-1,g,2);
assert(norm(pb(:)-b(:))/norm(b(:))<1e-10);
assert(norm(pr(:)-r(:))/norm(r(:))<1e-10);
fprintf('Independent exact measurement and physical Born/Rytov prediction: PASS\n');
end
