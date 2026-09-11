function low_contrast_xz
%LOW_CONTRAST_XZ Raw total Ex: exact Maxwell and scalar Rytov, five angles.
out=fileparts(mfilename('fullpath')); root=fileparts(fileparts(out));
addpath(fullfile(root,'forward','exp'),fullfile(root,'forward'));
api=oblate_na_sweep('helpers');
p=struct('lambda',.532,'n_m',1.335381534,'n_p',1.335381534+.01*(1.365-1.335381534), ...
    'a',3,'b',5,'NA_det',1,'z_det',5,'theta',[0;5;10;15;20]*pi/180);
p.illumination_NA=p.n_m*sin(p.theta);
x=linspace(-8,8,201); z=linspace(-5,9,176); [X,Z]=meshgrid(x,z);
points=[X(:),zeros(numel(X),1),Z(:)]; theta_deg=p.theta*180/pi;
k=2*pi*p.n_m/p.lambda;
u0=exp(1i*k*(points(:,1)*sin(p.theta.')+points(:,3)*cos(p.theta.')));
incidentEx=u0.*cos(p.theta.');
inside=(X(:)/p.b).^2+(Z(:)/p.a).^2<1-64*eps;
file=fullfile(out,'low_contrast_xz_fields.mat');
if isfile(file)
    old=load(file,'p'); assert(isequal(old.p,p),'XZ checkpoint parameter mismatch.');
    fprintf('LOW_CONTRAST_XZ_CACHED %s\n',file); return
end
s=api.reference(p,114,59);
[exactEx,scalar]=api.near_fields(s,points,true);
ids=unique(round(linspace(1,size(points,1),257)));
low=api.reference(p,109,55);
[el,ul]=api.near_fields(low,points(ids,:),true);
% Use the weak scattered-field norm, not the much larger incident-field norm.
vector_cutoff=relative(el-incidentEx(ids,:),exactEx(ids,:)-incidentEx(ids,:));
scalar_cutoff=relative(ul-u0(ids,:),scalar(ids,:)-u0(ids,:));
assert(max(vector_cutoff,scalar_cutoff)<1e-5,'Near-field cutoff convergence.');
clear low el ul

% Symmetric derivative at zero potential is the full-Green first Born field.
h=.005; plus=api.scalar_scaled_reference(s,h); minus=api.scalar_scaled_reference(s,-h);
derivative=plus;
for j=1:numel(s.modes)
    derivative.modes(j).so=(plus.modes(j).so-minus.modes(j).so)/(2*h);
end
born=complex(zeros(size(u0)));
[~,outside]=api.near_fields(derivative,points(~inside,:),false);
born(~inside,:)=outside-u0(~inside,:);
[~,up]=api.near_fields(plus,points(inside,:),false);
[~,um]=api.near_fields(minus,points(inside,:),false);
born(inside,:)=(up-um)/(2*h);
rytovEx=incidentEx.*exp(born./u0);

probe=[0,0,0;.7,0,.4;2,0,1;3.5,0,0;0,0,5;1,0,5;4,0,5;-6,0,2];
[~,up]=api.near_fields(plus,probe,false); [~,um]=api.near_fields(minus,probe,false);
bp=(up-um)/(2*h); clear plus minus derivative outside
plus=api.scalar_scaled_reference(s,2*h); minus=api.scalar_scaled_reference(s,-2*h);
[~,up]=api.near_fields(plus,probe,false); [~,um]=api.near_fields(minus,probe,false);
step_gap=relative((up-um)/(4*h),bp); clear plus minus
surface_gap=zeros(5,1);
for j=1:5
    [bs,~,~,info]=born_rytov_spheroid(p.lambda,p.n_p,p.n_m,p.a,p.b,p.theta(j),0, ...
        probe(:,1),probe(:,2),probe(:,3),struct('tol',1e-7));
    assert(info.validated); surface_gap(j)=relative(bp(:,j),bs);
end
assert(step_gap<1e-5 && max(surface_gap)<1e-5,'Born derivative checks.');
assert(all(isfinite(exactEx),'all') && all(isfinite(rytovEx),'all'));
validation=table(vector_cutoff,scalar_cutoff,step_gap,max(surface_gap), ...
    'VariableNames',{'vector_scattered_cutoff','scalar_scattered_cutoff', ...
    'born_step_gap','born_surface_gap'});
writetable(validation,fullfile(out,'low_contrast_xz_validation.csv'));
save(file,'p','x','z','theta_deg','u0','incidentEx','inside','exactEx','scalar','born','rytovEx','validation','-v7');
disp(validation);
fprintf('LOW_CONTRAST_XZ_PASS rows=%d angles=%d\n',size(points,1),numel(theta_deg));
end

function e=relative(a,b)
assert(norm(b(:))>0); e=norm(a(:)-b(:))/norm(b(:));
end
