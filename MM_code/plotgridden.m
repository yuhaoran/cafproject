SetDefault
clear; clc; close all
ng=128;
nplot=4;

fid=fopen('0.000den128_2.bin');
den=reshape(fread(fid,'real*4'),ng,ng,ng);
fclose(fid);

fid=fopen('0.000def_0128_2.bin');
def=reshape(fread(fid,'real*4'),ng,ng,ng);
fclose(fid);

% differential
dsp=zeros(3,ng,ng,ng);
for k=1:ng
for j=1:ng
for i=1:ng
  kp=mod(k,ng)+1;
  jp=mod(j,ng)+1;
  ip=mod(i,ng)+1;
  dsp(1,i,j,k)=def(ip,j,k)+def(ip,jp,k)+def(ip,j,kp)+def(ip,jp,kp)-def(i,j,k)-def(i,jp,k)-def(i,j,kp)-def(i,jp,kp);
  dsp(2,i,j,k)=def(i,jp,k)+def(ip,jp,k)+def(i,jp,kp)+def(ip,jp,kp)-def(i,j,k)-def(ip,j,k)-def(i,j,kp)-def(ip,j,kp);
  dsp(3,i,j,k)=def(i,j,kp)+def(ip,j,kp)+def(i,jp,kp)+def(ip,jp,kp)-def(i,j,k)-def(ip,j,k)-def(i,jp,k)-def(ip,jp,k);
end
end
end
dsp=dsp/4;

% get corrected den
k1=ng/2;
k2=k1+nplot-1;
den_corr=zeros(ng,ng);
for j=1:ng
for i=1:ng
  jm=mod(j-2,ng)+1;
  im=mod(i-2,ng)+1;
  k1m=mod(k1-2,ng)+1;
  %k2m=mod(k2-2,ng)+1;
  kk1=k1-1+(dsp(3,im,jm,k1m)+dsp(3,i,jm,k1m)+dsp(3,im,j,k1m)+dsp(3,i,j,k1m))/4;
  kk2=k2+(dsp(3,im,jm,k2)+dsp(3,i,jm,k2)+dsp(3,im,j,k2)+dsp(3,i,j,k2))/4;
  nk1=floor(kk1)+1;
  nk2=floor(kk2)+1;
  r1=nk1-kk1;
  r2=kk2-nk2+1;
  den_corr(i,j)=r1*den(i,j,nk1)+sum(den(i,j,nk1+1:nk2-1))+r2*den(i,j,nk2)-sum(den(i,j,nk2:nk1));
end
end
den_corr=den_corr/nplot;
%%
fig(1,'den');
imagesc(reshape(mean(den(:,:,k1:k2),3),ng,ng)')
colormap(hot)
caxis([0,1000*nplot])
axis xy square
colorbar
hold on

fid=fopen('newposition_128_2.dat');
newp=reshape(fread(fid,'real*4'),3,ng,ng,ng);
fclose(fid);
for jp=1:ng
    plot(reshape(mean(newp(1,:,jp,k1:k2),4),ng,1)',reshape(mean(newp(2,:,jp,k1:k2),4),ng,1)','Color',[0.5 0.5 0.5])
    axis xy square; hold on;
end
for ip=1:ng
    plot(reshape(mean(newp(1,ip,:,k1:k2),4),ng,1)',reshape(mean(newp(2,ip,:,k1:k2),4),ng,1)','Color',[0.5 0.5 0.5])
    hold on;
end
%%
fig(2,'den_corr');
imagesc(den_corr')
colormap(hot)
caxis([0,1000*nplot])
axis xy square
colorbar
hold on

fid=fopen('newposition_128_2.dat');
newp=reshape(fread(fid,'real*4'),3,ng,ng,ng);
fclose(fid);
for jp=1:ng
    plot(reshape(mean(newp(1,:,jp,k1:k2),4),ng,1)',reshape(mean(newp(2,:,jp,k1:k2),4),ng,1)','Color',[0.5 0.5 0.5])
    axis xy square; hold on;
end
for ip=1:ng
    plot(reshape(mean(newp(1,ip,:,k1:k2),4),ng,1)',reshape(mean(newp(2,ip,:,k1:k2),4),ng,1)','Color',[0.5 0.5 0.5])
    hold on;
end
