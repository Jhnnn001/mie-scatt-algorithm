function op = odt_operator(q, sz, dx, padding, weights)
%ODT_OPERATOR Matched padded FFT/trilinear Ewald operator and complex adjoint.
% op.forward maps chi to normalized weighted samples; physical Fourier
% integrals are op.forward(chi)./op.sample_scale. No q coordinates are rounded.
% Internal padding controls interpolation accuracy, not the direct FFT grid.
if nargin<4, padding=2; end
if nargin<5, weights=ones(size(q,1),1); end
validateattributes(q,{'numeric'},{'real','finite','nonempty','ncols',3});
validateattributes(sz,{'numeric'},{'integer','finite','numel',3,'>=',4});
validateattributes(dx,{'numeric'},{'real','finite','scalar','positive'});
validateattributes(padding,{'numeric'},{'integer','finite','scalar','positive'});
validateattributes(weights,{'numeric'},{'real','finite','vector','positive','numel',size(q,1)});
sz=double(sz(:).'); assert(all(mod(sz,2)==0),'Use even volume dimensions.');
ps=sz*padding; np=prod(ps); ns=size(q,1); rootn=sqrt(np);
weights=double(weights(:));
coords=q./(2*pi./(ps*dx))+ps/2+1;
coords=min(max(coords,1),ps);
bounds=[max(1,min(ps-1,floor(min(coords,[],1))));min(ps,floor(max(coords,[],1))+1)];
T=interpolation_transpose(q,ps,dx,false,bounds);
selected=arrayfun(@(j) bounds(1,j):bounds(2,j),1:3,'UniformOutput',false);
planes=selected{3}; spectral_size=bounds(2,:)-bounds(1,:)+1;
z=(-sz(3)/2:sz(3)/2-1)*dx;
kz=(planes-ps(3)/2-1)*2*pi/(ps(3)*dx);
axial=exp(-1i*z(:)*kz);
density=full(T*weights); bound=max(density);
rootweight=sqrt(weights/bound);
T=T*spdiags(rootweight,0,ns,ns);
% Direct gridding always uses the actual reconstruction grid. Its adjoint
% deposition is density normalized and zero filled, not an iterative fit.
S=interpolation_transpose(q,sz,dx,true);
S=S*spdiags(sqrt(weights),0,ns,ns);
direct_density=reshape(full(S*sqrt(weights)),sz);
direct_scale=sqrt(bound)*rootn/sqrt(prod(sz));
offset=(ps-sz)/2;
ix=offset(1)+(1:sz(1)); iy=offset(2)+(1:sz(2));
negative=arrayfun(@(s) mod(-(0:s-1),s)+1,sz,'UniformOutput',false);
op=struct('forward',@forward,'adjoint',@adjoint,'direct',@direct, ...
    'gridding',@gridding,'grid_size',sz,'padded_size',ps,'gridding_size',sz,'dx',dx, ...
    'padding',padding,'sample_scale',rootweight/(dx^3*rootn), ...
    'normal_diagonal',full(sum(sum(T.^2)))/np,'norm_bound',1, ...
    'physical_normalization',1/(bound*dx^6*np), ...
    'sample_count',ns,'interpolation','matched trilinear');
clear density q weights coords

    function y=forward(x)
        % Exact separable transform; pad each axis only when transforming it.
        padded=complex(zeros([ps(1),sz(2),numel(planes)]));
        padded(ix,:,:)=reshape(reshape(x,[],sz(3))*axial,[sz(1:2),numel(planes)]);
        along_x=fftshift(fft(ifftshift(padded,1),[],1),1);
        padded=complex(zeros([spectral_size(1),ps(2),numel(planes)]));
        padded(:,iy,:)=along_x(selected{1},:,:);
        ft=fftshift(fft(ifftshift(padded,2),[],2),2);
        ft=ft(:,selected{2},:)/rootn;
        y=T'*ft(:);
    end

    function x=adjoint(y)
        ft=reshape(T*y(:),spectral_size);
        padded=complex(zeros([spectral_size(1),ps(2),numel(planes)]));
        padded(:,selected{2},:)=ft;
        along_y=fftshift(ifft(ifftshift(padded,2),[],2),2)*ps(2);
        padded=complex(zeros([ps(1),sz(2),numel(planes)]));
        padded(selected{1},:,:)=along_y(:,iy,:);
        padded=fftshift(ifft(ifftshift(padded,1),[],1),1)*(ps(1)/rootn);
        x=reshape(reshape(padded(ix,:,:),[],numel(planes))*axial',sz);
    end

    function [ft,mask]=gridding(y)
        numerator=reshape(S*(y(:)*direct_scale),sz);
        denominator=direct_density;
        numerator=numerator+conj(numerator(negative{1},negative{2},negative{3}));
        denominator=denominator+denominator(negative{1},negative{2},negative{3});
        mask=denominator>0; ft=complex(zeros(sz));
        ft(mask)=numerator(mask)./denominator(mask);
    end

    function x=direct(y)
        ft=gridding(y);
        x=real(fftshift(ifftn(ifftshift(ft)))*sqrt(prod(sz)));
    end
end

function T=interpolation_transpose(q,ps,dx,periodic,bounds)
if nargin<5, bounds=[ones(1,3);ps]; end
spectral_size=bounds(2,:)-bounds(1,:)+1;
dk=2*pi./(ps*dx); coords=q./dk+ps/2+1; ns=size(q,1);
if periodic
    assert(all(abs(q)<pi/dx+1e-10,'all'),'Direct grid Nyquist violation.');
    coords=mod(coords-1,ps)+1; lo=floor(coords);
else
    if any(coords<1-1e-10 | coords>ps+1e-10,'all')
        error('odt_operator:Bandwidth','Ewald samples exceed the padded FFT bandwidth.');
    end
    coords=min(max(coords,1),ps); lo=min(floor(coords),ps-1);
end
fraction=coords-lo; columns=zeros(ns,8); values=columns;
for c=0:7
    bit=[bitget(c,1),bitget(c,2),bitget(c,3)];
    ix=mod(lo+bit-1,ps)+1;
    ix=ix-bounds(1,:)+1;
    columns(:,c+1)=ix(:,1)+(ix(:,2)-1)*spectral_size(1)+(ix(:,3)-1)*prod(spectral_size(1:2));
    values(:,c+1)=prod(bit.*fraction+(1-bit).*(1-fraction),2);
end
T=sparse(columns(:),repmat((1:ns).',8,1),values(:),prod(spectral_size),ns);
end
