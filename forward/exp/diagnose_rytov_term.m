function stats = diagnose_rytov_term(step)
% Evaluate the omitted complex-gradient square from an exact scalar solution.
% Q = grad(psi).grad(psi), NOT grad(psi)*conj(grad(psi)).
out=fileparts(mfilename('fullpath'));
if nargin==0, step=.1; end
addpath(fullfile(out,'..','..','spheroid-analytic-forward'));
lambda=.532; nm=1.335381534; np=1.37; a=5; b=2.5;
k=2*pi*nm/lambda; f0=(2*pi/lambda)^2*(np^2-nm^2);
rho=step/2:step:5-step/2; zz=-6+step/2:step:15-step/2;
[R,Z]=ndgrid(rho,zz); p=[R(:),zeros(numel(R),1),Z(:)];
inside=R.^2/b^2+Z.^2/a^2<1;
names={'interior_core','interior_rim','exit_5_to_7','downstream_7_to_14'};
masks={inside & R<b/2,inside & R>.8*b, ...
    Z>5 & Z<7 & R<3,Z>7 & Z<14 & R<3};
rows=struct([]); j=0;
for t=[.1,1]
    timer=tic;
    s=scalar_prolate(lambda,sqrt(nm^2+t*(np^2-nm^2)),nm,a,b,112);
    [u,du]=scalar_prolate_eval(s,p);
    dpsi=du./u; dpsi(:,3)=dpsi(:,3)-1i*k;
    Q=sum(dpsi.^2,2); transverse=dpsi(:,1).^2+dpsi(:,2).^2;
    finite=abs(u)>.1;
    % Independent finite-difference PDE check at off-boundary interior/exterior points.
    pp=[.3,0,0; 1.5,0,1; 2.1,0,2; 1,0,5.4; 2,0,7; 1,0,12];
    [v,dv]=scalar_prolate_eval(s,pp); lap=zeros(size(v)); h=2e-5;
    for dim=1:3
        plus=pp; minus=pp; plus(:,dim)=plus(:,dim)+h; minus(:,dim)=minus(:,dim)-h;
        [~,dp]=scalar_prolate_eval(s,plus); [~,dm]=scalar_prolate_eval(s,minus);
        lap=lap+(dp(:,dim)-dm(:,dim))/(2*h);
    end
    g=dv./v; g(:,3)=g(:,3)-1i*k; qp=sum(g.^2,2);
    lappsi=lap./v-sum((dv./v).^2,2);
    ff=t*f0*(sum(pp(:,1:2).^2,2)/b^2+pp(:,3).^2/a^2<1);
    residual=lappsi+2i*k*g(:,3)+qp+ff;
    pde_residual=norm(residual)/norm(ff);
    assert(pde_residual<1e-5,'Exact Rytov equation residual %.3g.',pde_residual);
    for region=1:numel(names)
        mask=masks{region}(:) & finite;
        w=sqrt(R(mask)); qabs=sort(abs(Q(mask)));
        j=j+1; rows(j).potential_scale=t; rows(j).region=string(names{region});
        rows(j).Q_rms_um_minus2=norm(w.*Q(mask))/norm(w);
        rows(j).Q_median_um_minus2=qabs(ceil(end*.5));
        rows(j).Q_p90_um_minus2=qabs(ceil(end*.9));
        rows(j).Q_over_f_rms=NaN;
        if region<=2, rows(j).Q_over_f_rms=rows(j).Q_rms_um_minus2/(t*f0); end
        rows(j).transverse_fraction=norm(w.*transverse(mask))/norm(w.*Q(mask));
        rows(j).min_abs_U=min(abs(u(masks{region}(:))));
        rows(j).masked_fraction=1-nnz(mask)/nnz(masks{region});
        rows(j).equation_residual=pde_residual;
        fprintf('t=%.1f %-20s rms|Q|=%.4g, Q/f=%.4g, transverse ratio=%.4g, min|U|=%.4g\n', ...
            t,names{region},rows(j).Q_rms_um_minus2,rows(j).Q_over_f_rms, ...
            rows(j).transverse_fraction,rows(j).min_abs_U);
    end
    save(fullfile(out,sprintf('term_t%02d_dx%g.mat',round(t*10),step)), ...
        'R','Z','inside','u','dpsi','Q','transverse','finite','f0','t','pde_residual');
    fprintf('Term map t=%.1f: PDE residual %.3g, %.2fs.\n',t,pde_residual,toc(timer));
end
stats=struct2table(rows);
writetable(stats,fullfile(out,sprintf('omitted_term_dx%g.csv',step)));
% Infinite slab control: identical index contrast and 10 um central thickness.
kp=2*pi*np/lambda; d=2*a; r=(k-kp)/(k+kp);
T=(1-r^2)*exp(1i*(kp-k)*d)/(1-r^2*exp(2i*kp*d));
R=exp(1i*f0*d/(2*k));
slab_total_error=abs(R-T)/abs(T); slab_scattered_error=abs(R-T)/abs(T-1);
fprintf('Uniform 10 um slab: total error %.4f%%, scattered error %.4f%%, phase %.4f rad.\n', ...
    100*slab_total_error,100*slab_scattered_error,(kp-k)*d);
save(fullfile(out,'slab_control.mat'),'T','R','slab_total_error','slab_scattered_error');
end
