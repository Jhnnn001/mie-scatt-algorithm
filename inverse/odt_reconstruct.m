function [n, chi, info] = odt_reconstruct(data, method, opts)
%ODT_RECONSTRUCT Direct Fourier, nonnegative GP, or Lim TV reconstruction.
% data.operator is an odt_operator and data.y its normalized weighted samples.
% TV: min alpha*TV(chi)+.5*||A chi-g_k||^2, chi>=0, followed by
% g_{k+1}=g_k+g-A chi (Lim 2015 Eqs.14-15 and Appendix A).
% outer_iter=1 is fixed-data penalized TV. Differences use unit voxel pitch.
% GP uses a CGLS minimum-norm data correction, then positivity; projection
% accuracy is reported separately. Inconsistent data use least-squares fits.
started=tic;
if nargin<3, opts=struct; end
method=validatestring(method,{'direct','gp','tv'});
defaults=struct('alpha',1e-4,'rho',.1,'rho_positive',.1,'max_iter',100, ...
    'outer_iter',5,'tol',1e-4,'cg_tol',1e-3,'cg_max_iter',50, ...
    'discrepancy',0,'verbose',false);
if ~isstruct(opts) || ~isscalar(opts) || ...
        any(~ismember(fieldnames(opts),fieldnames(defaults)))
    error('odt_reconstruct:InvalidOptions','Unknown solver options.');
end
for name=fieldnames(defaults).'
    if ~isfield(opts,name{1}), opts.(name{1})=defaults.(name{1}); end
end
for name={'alpha','discrepancy'}
    validateattributes(opts.(name{1}),{'numeric'},{'real','finite','scalar','nonnegative'});
end
for name={'rho','rho_positive','tol','cg_tol'}
    validateattributes(opts.(name{1}),{'numeric'},{'real','finite','scalar','positive'});
end
for name={'max_iter','outer_iter','cg_max_iter'}
    validateattributes(opts.(name{1}),{'numeric'},{'integer','finite','scalar','positive'});
end
validateattributes(opts.verbose,{'logical'},{'scalar'});
if ~isstruct(data) || ~all(isfield(data,{'operator','y','n_m'}))
    error('odt_reconstruct:InvalidData','Expected operator, y and n_m.');
end
validateattributes(data.n_m,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(data.y,{'numeric'},{'finite','vector'});
op=data.operator; sz=op.grid_size; y=double(data.y(:));
assert(numel(y)==op.sample_count,'Data size does not match operator.');
A=op.forward; AH=@(r) real(op.adjoint(r));
chi=op.direct(y);
alpha=opts.alpha*strcmp(method,'tv');
initial=chi; if ~strcmp(method,'direct'), initial=max(initial,0); end
[initial_objective,~]=objective(initial,y,alpha,A);
history=zeros(opts.max_iter*opts.outer_iter,8);
outer_history=zeros(opts.outer_iter,5);
iteration=0; outer=0; reinjection_norm=0;
converged=strcmp(method,'direct'); inner_converged=converged;
projection_converged=true; cg_failures=0;
stationarity=NaN;
ynorm=max(norm(y),eps); aty=AH(y);
normal_scale=max(norm(aty(:)),eps);
if norm(y)==0
    chi=zeros(sz); converged=true; inner_converged=true;
    stationarity=0;
elseif strcmp(method,'gp')
    chi=max(chi,0); outer=1;
    for k=1:opts.max_iter
        old=chi;
        [correction,projection_res,cgit]=cgls(y-A(chi),A,AH,sz, ...
            opts.cg_tol*normal_scale,opts.cg_max_iter);
        projected=chi+correction;
        chi=max(projected,0);
        change=norm(chi(:)-old(:))/max(norm(old(:)),eps);
        [cost,residual]=objective(chi,y,0,A);
        iteration=iteration+1;
        projection_converged=projection_res<=opts.cg_tol*normal_scale;
        cg_failures=cg_failures+~projection_converged;
        history(iteration,:)=[outer,k,cost,residual,change,0,projection_res/normal_scale,cgit];
        if opts.verbose && (mod(k,10)==0 || k==1)
            fprintf('GP %d: data=%.5g step=%.4g projection=%.4g CG=%d time=%.1fs\n', ...
                k,residual,change,projection_res/normal_scale,cgit,toc(started));
        end
        if change<=opts.tol && projection_converged
            converged=true; break;
        end
    end
    inner_converged=converged;
elseif strcmp(method,'tv')
    chi=max(chi,0); v=chi; bv=zeros(sz);
    d=difference(chi); b=zeros([sz,3]); gk=y;
    lap=zeros(sz);
    for axis=1:3
        dims=ones(1,3); dims(axis)=sz(axis);
        lap=lap+reshape(4*sin(pi*(-sz(axis)/2:sz(axis)/2-1)/sz(axis)).^2,dims);
    end
    denominator=op.normal_diagonal+opts.rho*lap+opts.rho_positive;
    precondition=@(r) preconditioner(r,denominator,sz);
    normal=@(x) normal_equation(x,A,AH,opts,sz);
    for outer=1:opts.outer_iter
        atg=AH(gk); inner_converged=false;
        for k=1:opts.max_iter
            rhs=atg+opts.rho*adjoint_difference(d-b)+opts.rho_positive*(v-bv);
            % Scale the solve to the data, not a possibly huge penalty RHS.
            cg_tolerance=min(opts.cg_tol,.1*opts.tol*normal_scale/max(norm(rhs(:)),eps));
            [solution,cgflag,cgres,cgit]=pcg(normal,rhs(:),cg_tolerance, ...
                opts.cg_max_iter,precondition,[],chi(:));
            cg_failures=cg_failures+(cgflag~=0);
            chi=reshape(solution,sz); grad=difference(chi);
            old_d=d; old_v=v;
            z=grad+b; magnitude=sqrt(sum(z.^2,4));
            d=z.*max(1-(opts.alpha/opts.rho)./max(magnitude,realmin),0);
            v=max(chi+bv,0); b=b+grad-d; bv=bv+chi-v;
            primal=hypot(norm(grad(:)-d(:)),norm(chi(:)-v(:)))/ ...
                max([hypot(norm(grad(:)),norm(chi(:))),hypot(norm(d(:)),norm(v(:))),eps]);
            dual_field=opts.rho*adjoint_difference(d-old_d)+opts.rho_positive*(v-old_v);
            dual=norm(dual_field(:))/normal_scale;
            [cost,residual]=objective(v,y,opts.alpha,A);
            iteration=iteration+1;
            history(iteration,:)=[outer,k,cost,residual,primal,dual,cgres,cgit];
            if opts.verbose && mod(k,100)==0
                fprintf('TV %d/%d: data=%.5g primal=%.4g dual=%.4g CG=%d time=%.1fs\n', ...
                    outer,k,residual,primal,dual,cgit,toc(started));
            end
            if max(primal,dual)<=opts.tol && cgflag==0
                kkt=AH(A(v)-gk)+opts.rho*adjoint_difference(b)+opts.rho_positive*bv;
                stationarity=norm(kkt(:))/normal_scale;
                if stationarity<=opts.tol
                    inner_converged=true; break;
                end
            end
        end
        chi=v;
        kkt=AH(A(chi)-gk)+opts.rho*adjoint_difference(b)+opts.rho_positive*bv;
        stationarity=norm(kkt(:))/normal_scale;
        residual=norm(A(chi)-y)/ynorm;
        outer_history(outer,:)=[outer,k,residual,primal,dual];
        if opts.verbose
            fprintf('TV outer %d inner %d: data=%.5g primal=%.4g dual=%.4g time=%.1fs\n', ...
                outer,k,residual,primal,dual,toc(started));
        end
        target=max(opts.discrepancy,opts.tol);
        if inner_converged && (opts.outer_iter==1 || residual<=target)
            converged=true; break;
        end
        if outer<opts.outer_iter
            increment=y-A(chi);
            gk=gk+increment;
            reinjection_norm=reinjection_norm+norm(increment);
        end
    end
end
if any(~isfinite(chi(:))) || any(chi(:)<-1)
    error('odt_reconstruct:InvalidResult','Nonfinite chi or negative n squared.');
end
n=data.n_m*sqrt(1+chi);
[cost,residual]=objective(chi,y,alpha,A);
if iteration==0
    history=[0,0,cost,residual,0,0,0,0];
else
    history=history(1:iteration,:);
end
info=struct('method',method,'iterations',iteration,'outer_iterations',outer, ...
    'converged',converged,'inner_converged',inner_converged, ...
    'projection_converged',projection_converged,'cg_failures',cg_failures, ...
    'projection_criterion','CGLS normal residual, not a solution-error bound', ...
    'stationarity_residual',stationarity, ...
    'data_consistent',residual<=opts.tol,'relative_data_residual',residual, ...
    'objective',history(:,3),'residual',history(:,4), ...
    'primal_or_step_residual',history(:,5),'dual_residual',history(:,6), ...
    'history',history,'outer_history',outer_history(1:outer,:), ...
    'initial_objective',initial_objective,'alpha',alpha,'options',opts, ...
    'reinjection_norm',reinjection_norm,'time',toc(started),'boundary','periodic');
end

function [x,residual,iteration]=cgls(r,A,AH,sz,tolerance,maxiter)
% Zero-start CGLS gives the minimum-norm correction, to reported tolerance.
x=zeros(sz); s=AH(r); p=s; gamma=sum(s(:).^2);
residual=sqrt(gamma); iteration=0;
for k=1:maxiter
    if residual<=tolerance, break; end
    q=A(p); denominator=sum(abs(q).^2);
    if denominator<=realmin, break; end
    step=gamma/denominator; x=x+step*p; r=r-step*q;
    s=AH(r); next=sum(s(:).^2);
    p=s+(next/gamma)*p; gamma=next; residual=sqrt(gamma); iteration=k;
end
end

function y=normal_equation(x,A,AH,opts,sz)
x=reshape(x,sz);
y=AH(A(x))+opts.rho*adjoint_difference(difference(x))+opts.rho_positive*x;
y=y(:);
end

function x=preconditioner(r,denominator,sz)
x=real(fftshift(ifftn(ifftshift(fftshift(fftn(ifftshift(reshape(r,sz))))./denominator))));
x=x(:);
end

function grad=difference(x)
grad=zeros([size(x),3]);
for axis=1:3, grad(:,:,:,axis)=circshift(x,-1,axis)-x; end
end

function x=adjoint_difference(grad)
x=zeros(size(grad,1),size(grad,2),size(grad,3));
for axis=1:3
    part=grad(:,:,:,axis); x=x+circshift(part,1,axis)-part;
end
end

function [cost,residual]=objective(x,y,alpha,A)
r=A(x)-y; cost=.5*sum(abs(r).^2);
if alpha>0
    grad=difference(x); cost=cost+alpha*sum(sqrt(sum(grad.^2,4)),'all');
end
residual=norm(r)/max(norm(y),eps);
end
