function root = locate_field_zero
% A local zero of the pre-objective total field makes log(U/U0) singular.
out=fileparts(mfilename('fullpath'));
addpath(fullfile(out,'..','..','spheroid-analytic-forward'));
m=load(fullfile(out,'term_t10_dx0.025.mat'));
[~,ix]=min(abs(m.u)); p=[m.R(ix),0,m.Z(ix)];
s=scalar_prolate(.532,1.37,1.335381534,5,2.5,112);
for iteration=1:15
    [u,du]=scalar_prolate_eval(s,p);
    if abs(u)<1e-11, break; end
    J=[real(du([1,3]));imag(du([1,3]))];
    shift=J\[real(u);imag(u)]; p([1,3])=p([1,3])-shift.';
end
assert(p(1)>0 && p(3)>5 && abs(u)<1e-9,'Local field-zero solve failed.');
refined=scalar_prolate(.532,1.37,1.335381534,5,2.5,132);
check=scalar_prolate_eval(refined,p);
assert(abs(check)<1e-8,'Field zero did not survive L refinement.');
root=struct('rho_um',p(1),'z_um',p(3),'abs_U',abs(u),'abs_U_refined',abs(check));
writetable(struct2table(root),fullfile(out,'field_zero.csv'));
fprintf('Pre-objective field zero: rho=%.10g z=%.10g um; |U|=%.3g, refined %.3g.\n', ...
    p(1),p(3),abs(u),abs(check));
end
