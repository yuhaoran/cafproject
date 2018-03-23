clc; clear;
SetDefault
nn=1;
nnt=2;
nc=64;
ng=nc*4;
nfb=24;
nt=nc/nnt;
nft=4*nt;
nfe=nft+2*nfb;
npen=nc/nn;
ng_global=ng*nn;
nplocal=ng_global^3/8;
box=200;
Pk_shot=box^3/nplocal;
universe=2;
%close all
%% HMF
ninfo=20;
fid=fopen(['../../cafproject/CUBEnu/output/universe',int2str(universe),'/image1/0.000_halo_1.bin']);
nhalo_tot=fread(fid,1,'integer*4')';
nhalo=fread(fid,1,'integer*4')';
halo_vir=fread(fid,1,'real*4')';
halo_odc=fread(fid,1,'real*4')';
haloinfo=fread(fid,[ninfo,nhalo],'real*4');

%% quick check power
figure

n_row_xi=10;
fid=fopen(['../../cafproject/CUBEnu/output/universe',int2str(universe),'/image1/0.000_cicpower_1.bin']);
xi=fread(fid,'real*4')';
fclose(fid);
xi=reshape(xi,n_row_xi,numel(xi)/n_row_xi)';

subplot(2,2,1)
a=load('CAMB/pkz0.txt');
%loglog(a(:,1),a(:,2).*a(:,1).^-0.04*0.88,'--'); hold on
k3=xi(:,2).^3/(2*pi^2);
Deltak_shot=Pk_shot*k3;
%loglog(xi(:,2),xi(:,3)./k3-Pk_shot,xi(:,2),xi(:,4)./k3-Pk_shot,xi(:,2),xi(:,5)./k3-Pk_shot)
loglog(a(:,1),a(:,2).*a(:,1).^3/(2*pi^2),'--',...
       xi(:,2),(xi(:,3)-Deltak_shot),...
       xi(:,2),Deltak_shot,'--')
grid on; hold on
label('$k$/($h$ Mpc$^{-1}$)','$\Delta^2(k)$')
legend('CLASS','$cc-$shot','shot','Location','northwest')
ax1=gca;
  ax1.XLim=[xi(1,2),xi(length(xi),2)];
  ax1.YLim=[1e-3,1e3];

subplot(2,2,3)
semilogx(xi(:,2),xi(:,8))
grid on; hold on
label('$k$/($h$ Mpc$^{-1}$)','$\xi_{NL}$')
ax1=gca;
  ax1.XLim=[xi(1,2),xi(length(xi),2)];
  ax1.YLim=[0,1];
  
subplot(2,2,2)
a=loadfield3d(['../../cafproject/CUBEnu/output/universe',int2str(universe),'/image1/0.000_delta_c_1.bin']);
imagesc(reshape(mean(a(:,:,:),3),ng,ng)'); hold on
axis xy square; colorbar; colormap(1-gray); caxis([-1,5]); title('$\delta_c$')
plot(haloinfo(1,:),haloinfo(2,:),'r.')

subplot(2,2,4)
a=loadfield3d(['../../cafproject/CUBEnu/output/universe',int2str(universe),'/image1/delta_L_1.bin']);
imagesc(reshape(mean(a(:,:,:),3),ng,ng)');
axis xy square; colorbar; colormap(1-gray)
caxis([-1,5])
title('$\delta_L$')
return

%% HMF
ninfo=20;
fid=fopen(['../../cafproject/CUBEnu/output/universe',int2str(universe),'/image1/0.000_halo_1.bin']);
nhalo_tot=fread(fid,1,'integer*4')';
nhalo=fread(fid,1,'integer*4')';
halo_vir=fread(fid,1,'real*4')';
halo_odc=fread(fid,1,'real*4')';

haloinfo=fread(fid,[ninfo,nhalo],'real*4');
figure
plot(haloinfo(1,:),haloinfo(2,:),'.')

[hmass_vir,hcount_vir]=hist(haloinfo(4,:));
[hmass_odc,hcount_odc]=hist(haloinfo(5,:));
figure
loglog(hcount_vir,hmass_vir,hcount_odc,hmass_odc)
legend('vir','odc')
fclose(fid);
