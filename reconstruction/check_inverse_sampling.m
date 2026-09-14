function check_inverse_sampling
%CHECK_INVERSE_SAMPLING Separate continuous-model, voxel and interpolation errors.
here=fileparts(mfilename('fullpath')); addpath(here); out=fullfile(here,'results');
loaded=load(fullfile(out,'exact_detector_512.mat'),'raw'); raw=loaded.raw;
nm=raw.parameters.n_m; np=raw.parameters.n_p; chi0=np^2/nm^2-1;
q=raw.q; Q=sqrt(25*(q(:,1).^2+q(:,2).^2)+9*q(:,3).^2);
form=ones(size(Q)); nz=abs(Q)>1e-5;
form(nz)=3*(sin(Q(nz))-Q(nz).*cos(Q(nz)))./Q(nz).^3;
continuous=chi0*100*pi*form;
weights=sqrt(repelem(raw.weights,nnz(raw.geometry.pupil)));
physical=raw.samples;
continuous_error=norm((continuous-physical).*weights)/norm(physical.*weights);
fprintf('Continuous true-object Ewald mismatch %.7f%%\n',100*continuous_error);
configs={2,.1,[128,128,80];3,.1,[128,128,80];4,.1,[128,128,80]; ...
    6,.1,[128,128,80];8,.1,[128,128,80];3,.08,[160,160,100];6,.08,[160,160,100]};
rows=struct([]);
for j=1:size(configs,1)
    [padding,dx,sz]=configs{j,:};
    [X,Y,Z]=ndgrid((-sz(1)/2:sz(1)/2-1)*dx, ...
        (-sz(2)/2:sz(2)/2-1)*dx,(-sz(3)/2:sz(3)/2-1)*dx);
    chi=chi0*double((X/5).^2+(Y/5).^2+(Z/3).^2<1); clear X Y Z
    timer=tic; op=odt_operator(q,sz,dx,padding,weights.^2);
    predicted=op.forward(chi)./op.sample_scale;
    row=struct('padding',padding,'dx',dx,'nx',sz(1),'nz',sz(3), ...
        'continuous_to_physical',continuous_error, ...
        'voxel_to_continuous',norm((predicted-continuous).*weights)/norm(continuous.*weights), ...
        'voxel_to_physical',norm((predicted-physical).*weights)/norm(physical.*weights), ...
        'seconds',toc(timer));
    rows=[rows;row]; %#ok<AGROW>
    writetable(struct2table(rows),fullfile(out,'operator_convergence.csv'));
    fprintf('Padding %d dx=%g: voxel/continuous %.6f%% physical %.6f%% %.1fs\n', ...
        padding,dx,100*row.voxel_to_continuous,100*row.voxel_to_physical,row.seconds);
    clear op chi predicted
end
% Change detector pitch/window independently, retaining the exact same q's.
indices=[1,6,11,16,21,26,31]; ns=nnz(raw.geometry.pupil);
baseline=reshape(raw.samples,ns,81); baseline=baseline(:,indices);
[BX,BY]=ndgrid(raw.geometry.kx,raw.geometry.ky);
bx=BX(raw.geometry.pupil); by=BY(raw.geometry.pupil);
rows=struct([]);
for dx=[.05,.1]
    m=odt_exact_data(1024,dx,indices);
    d=odt_prepare(m.U,m.u0,m.geometry,'rytov',struct('weights',m.weights,'padding',1));
    g=m.geometry; dk=g.kx(2)-g.kx(1);
    ix=round((bx-g.kx(1))/dk)+1; iy=round((by-g.ky(1))/dk)+1;
    lookup=zeros(size(g.pupil)); lookup(g.pupil)=1:nnz(g.pupil);
    selected=lookup(sub2ind(size(g.pupil),ix,iy)); assert(all(selected>0));
    fine=reshape(d.samples,nnz(g.pupil),numel(indices)); fine=fine(selected,:);
    w=sqrt(raw.weights(indices).');
    numerator=(fine-baseline).*w; denominator=baseline.*w;
    row=struct('detector_n',1024,'detector_dx',dx,'window_um',1024*dx, ...
        'angles',numel(indices),'rytov_data_relative_gap',norm(numerator(:))/norm(denominator(:)));
    rows=[rows;row]; %#ok<AGROW>
    writetable(struct2table(rows),fullfile(out,'detector_convergence.csv'));
    fprintf('Detector dx=%g window=%g: Rytov data gap %.7f%%\n',dx,1024*dx,100*row.rytov_data_relative_gap);
    clear m d fine
end
fprintf('INVERSE_SAMPLING_CHECKS_PASS\n');
end
