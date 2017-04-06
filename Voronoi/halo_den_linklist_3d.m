%% 
SetDefault; clear; clc
ng=256;
load xhalo.mat
nhalo=length(xhalo);

% rescale
xhalo=xhalo*ng/4608;
if 0 % randomize halo positions
  xhalo=rand(3,nhalo)*ng;
end
area_halo=ones(1,nhalo);
den=zeros(ng,ng,ng);
den_ngp=den;
box=ng;
mesh_phy=[box/ng/2,box-box/ng/2];
%% 1. hoc_grid, ll_p
tic
hoc_g=zeros(ng,ng,ng);
ll_p=zeros(1,nhalo);
for ihalo=1:nhalo
  i=floor(xhalo(1,ihalo))+1;
  j=floor(xhalo(2,ihalo))+1;
  k=floor(xhalo(3,ihalo))+1;
  den_ngp(i,j,k)=den_ngp(i,j,k)+1;
  ll_p(ihalo)=hoc_g(i,j,k);
  hoc_g(i,j,k)=ihalo;
end
%plot(xhalo(1,:),xhalo(2,:),'r.'); grid on
%sum(sum((hoc_g==0)))
%imagesc(mean(den_ngp(:,:,:),3))
%return
%% 2. hoc_p, ll_grid
hoc_p=zeros(1,nhalo);
ll_g=zeros(1,ng^3);
for k=1:ng
for j=1:ng
for i=1:ng
  ig=ng^2*(k-1)+ng*(j-1)+(i-1)+1;
  % loop over non-empty grids only
  pp=hoc_g(i,j,k);
  while pp~=0
    %ll_g(ig)=hoc_p(pp); % assigns zero actually, don't need this step
    hoc_p(pp)=ig;
    pp=ll_p(pp);
  end
end
end
end

for k=1:ng
disp(k)
for j=1:ng
for i=1:ng
  ig=ng^2*(k-1)+ng*(j-1)+(i-1)+1;
  % loop over empty grids only
  if hoc_g(i,j,k)==0 % find nearest particle
    r=0; % initialize search radius
    search=1; % if we need to search neighboors
    r2min=1000*ng^2;
    while search==1 % search for (1+2r)^3 grids
      gpos=[i-0.5;j-0.5;k-0.5];
      r=r+1;
      for k_s=k-r:k+r
      for j_s=j-r:j+r
      for i_s=i-r:i+r
        ii_s=mod(i_s-1,ng)+1;
        jj_s=mod(j_s-1,ng)+1;
        kk_s=mod(k_s-1,ng)+1;
        pp=hoc_g(ii_s,jj_s,kk_s);
        while pp~=0 % found particle in current grid
          search=0; % cancel next searching radius
          hpos=xhalo(:,pp);
          dpos=hpos-gpos;
          dpos=mod(dpos+ng/2,ng)-ng/2;
          r2=sum(dpos.^2);
          if r2<r2min % check if the particle is nearest
            r2min=r2;
            %hoc_extended(i,j)=pp; % assign the nearest particle to current grid
            pp_candidate=pp;
          end
          pp=ll_p(pp);
        end % finished searching current grid
      end
      end
      end
    end % finished searching current layer
    % found nearest particle, exit searching loop
    ll_g(ig)=hoc_p(pp_candidate); % let ll_g point to pp_candidate's chain
    hoc_p(pp_candidate)=ig; % replace pp_candidate's hoc to be current grid
    area_halo(pp_candidate)=area_halo(pp_candidate)+1;
  end
end
end
end
%sum(sum(ll_g(1,:,:)~=0))

%% 3. assign density
for ihalo=1:nhalo
  ig=hoc_p(ihalo);
  while ig~=0
    % assign
    idx(3)=floor((ig-1)/ng^2)+1;
    idx(2)=floor((ig-1-ng^2*(idx(3)-1))/ng)+1;
    idx(1)=mod(ig-1,ng)+1;
    if (ig~=(idx(3)-1)*ng^2+(idx(2)-1)*ng+idx(1))
      stop
    end
    den(idx(1),idx(2),idx(3))=den(idx(1),idx(2),idx(3))+1/area_halo(ihalo);
    ig=ll_g(ig);
  end
end
toc
disp('sum,min')
disp(sum(sum(sum(den))))
disp(min(min(min(den))))
sum(sum(sum(den==0)))
