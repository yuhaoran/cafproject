function [ sim ] = zip2xv( prefix )
%% read zip2
fid=fopen([prefix,'zip2_0.dat']);
sim.nplocal=fread(fid,1,'integer*4');
sim.a=fread(fid,1,'real*4');
sim.t=fread(fid,1,'real*4');
sim.tau=fread(fid,1,'real*4');
sim.nts=fread(fid,1,'integer*4');

sim.dt=fread(fid,3,'integer*4');
sim.nstep=fread(fid,3,'integer*4');

sim.mass_p=fread(fid,1,'real*4');
sim.v_i2r=1./fread(fid,3,'real*4');
sim.shake_offset=fread(fid,3,'real*4');

sim.box=fread(fid,1,'real*4');
sim.rank=fread(fid,1,'integer*4');
sim.nn=fread(fid,1,'integer*2');
sim.nnt=fread(fid,1,'integer*2');
sim.nt=fread(fid,1,'integer*2');
sim.ncell=fread(fid,1,'integer*2');
sim.ncb=fread(fid,1,'integer*2');
sim.izipx=fread(fid,1,'integer*1');
sim.izipv=fread(fid,1,'integer*1');

sim.h0=fread(fid,1,'real*4');
sim.omega_m=fread(fid,1,'real*4');
sim.omega_l=fread(fid,1,'real*4');
sim.s8=fread(fid,1,'real*4');

sim.m_neu=fread(fid,3,'real*4');
sim.vsim2phys=fread(fid,1,'real*4');
sim.z_i=fread(fid,1,'real*4');

nc=sim.nt*sim.nnt;
nf=nc*sim.ncell;
%rhoc=eval(['fread(fid,nc^3,''integer*',int2str(sim.izip2),''')']);
rhoc=fread(fid,nc^3,'integer*4');
fclose(fid);

%if sim.izip2==1
%    fid3=fopen([prefix,'zip3_0.dat']);
%    temp=fread(fid3,'integer*4');
%    close(fid3)
%    rhoc(rhoc==255)=temp;
%end
%clear temp

%if sim.tile_struc==0
%    rhoc=reshape(rhoc,nc,nc,nc);
%else
    rhoc=reshape(rhoc,sim.nt,sim.nt,sim.nt,sim.nnt,sim.nnt,sim.nnt);
%end

%% read zip0,1
fid=fopen([prefix,'zip0_0.dat']);
xzip=eval(['fread(fid,[3,sim.nplocal],''uint',int2str(8*sim.izipx),''')']);
fclose(fid);
%fid=fopen([prefix,'zip1_0.dat']);
%vzip=eval(['fread(fid,[3,sim.nplocal],''integer*',int2str(sim.izipv),''')']);
%fclose(fid);
%% xv
sim.xv=zeros(6,sim.nplocal);
ip=0;
%if sim.tile_struc==0
%    for k=1:nc
%    for j=1:nc
%    for i=1:nc
%        np=rhoc(i,j,k);
%        for l=1:np
%            ip=ip+1;
%            sim.xv(1:3,ip)=([i;j;k]-1)*ncell + (xzip(:,ip)+0.5)/2^(izipx*8)*ncell;
%        end
%   end
%    end
%    end
%else
    for itz=1:sim.nnt
    for ity=1:sim.nnt
    for itx=1:sim.nnt
        for k=1:sim.nt
        for j=1:sim.nt
        for i=1:sim.nt
            np=rhoc(i,j,k,itx,ity,itz);
            for l=1:np
                ip=ip+1;
                sim.xv(1:3,ip)=([itx;ity;itz]-1)*sim.nt*4 + ([i;j;k]-1)*4 + (xzip(:,ip)+0.5)/2^(sim.izipx*8)*4;
            end
        end
        end
        end
    end
    end
    end
%end
if 1
    figure
    plot3(sim.xv(1,:),sim.xv(2,:),sim.xv(3,:),'k.'); axis square xy
    ax1=gca;
    ax1.XLim=[0,nf];
    ax1.YLim=[0,nf];
    ax1.ZLim=[0,nf];
end
%% read zipid
fid=fopen([prefix,'zipid_0.dat']);
sim.iczip=eval(['fread(fid,[4,sim.nplocal],''integer*',int2str(2),''')']);
fclose(fid);

sim.ic=nf*(sim.iczip(2:4,:)+32768.5)/65536;
sim.disp=mod(sim.xv(1:3,:)-sim.ic+nf/2,nf)-nf/2;
return
if 1
    figure
    plot3(sim.disp(1,:),sim.disp(2,:),sim.disp(3,:),'k.')
    grid on
end
%% density
n=nf*2;
sim.rho_f=zeros(n+2,n+2,n+2);
mass_p=sim.mass_p;
for ip=1:sim.nplocal
    tempx=sim.xv(1:3,ip)*n/nf+0.5;
    idx1=floor(tempx)+1;
    idx2=idx1+1;
    dx1=idx1-tempx;
    dx2=1-dx1;
    sim.rho_f(idx1(1),idx1(2),idx1(3))=sim.rho_f(idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*mass_p;
    sim.rho_f(idx2(1),idx1(2),idx1(3))=sim.rho_f(idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*mass_p;
    sim.rho_f(idx1(1),idx2(2),idx1(3))=sim.rho_f(idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*mass_p;
    sim.rho_f(idx1(1),idx1(2),idx2(3))=sim.rho_f(idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*mass_p;
    sim.rho_f(idx1(1),idx2(2),idx2(3))=sim.rho_f(idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*mass_p;
    sim.rho_f(idx2(1),idx1(2),idx2(3))=sim.rho_f(idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*mass_p;
    sim.rho_f(idx2(1),idx2(2),idx1(3))=sim.rho_f(idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*mass_p;
    sim.rho_f(idx2(1),idx2(2),idx2(3))=sim.rho_f(idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*mass_p;
end
sim.rho_f(2,:,:)=sim.rho_f(2,:,:)+sim.rho_f(n+2,:,:);
sim.rho_f(n+1,:,:)=sim.rho_f(n+1,:,:)+sim.rho_f(1,:,:);
sim.rho_f(:,2,:)=sim.rho_f(:,2,:)+sim.rho_f(:,n+2,:);
sim.rho_f(:,n+1,:)=sim.rho_f(:,n+1,:)+sim.rho_f(:,1,:);
sim.rho_f(:,:,2)=sim.rho_f(:,:,2)+sim.rho_f(:,:,n+2);
sim.rho_f(:,:,n+1)=sim.rho_f(:,:,n+1)+sim.rho_f(:,:,1);

sim.rho_f=sim.rho_f(2:n+1,2:n+1,2:n+1);
if 1
    figure;
    imagesc(reshape(sum(sim.rho_f,3),n,n)'); axis xy square; colorbar
end
%% disp field
n=nf/2;
sim.disp=zeros(3,n+2,n+2,n+2);
for ip=1:sim.nplocal
    tempx=sim.ic(1:3,ip)*n/nf+0.5;
    idx1=floor(tempx)+1;
    idx2=idx1+1;
    dx1=idx1-tempx;
    dx2=1-dx1;
    dr=mod(sim.xv(1:3,ip)-sim.ic(1:3,ip)+nf/2,nf)-nf/2;
    sim.disp(:,idx1(1),idx1(2),idx1(3))=sim.disp(:,idx1(1),idx1(2),idx1(3))+dx1(1)*dx1(2)*dx1(3)*dr;
    sim.disp(:,idx2(1),idx1(2),idx1(3))=sim.disp(:,idx2(1),idx1(2),idx1(3))+dx2(1)*dx1(2)*dx1(3)*dr;
    sim.disp(:,idx1(1),idx2(2),idx1(3))=sim.disp(:,idx1(1),idx2(2),idx1(3))+dx1(1)*dx2(2)*dx1(3)*dr;
    sim.disp(:,idx1(1),idx1(2),idx2(3))=sim.disp(:,idx1(1),idx1(2),idx2(3))+dx1(1)*dx1(2)*dx2(3)*dr;
    sim.disp(:,idx1(1),idx2(2),idx2(3))=sim.disp(:,idx1(1),idx2(2),idx2(3))+dx1(1)*dx2(2)*dx2(3)*dr;
    sim.disp(:,idx2(1),idx1(2),idx2(3))=sim.disp(:,idx2(1),idx1(2),idx2(3))+dx2(1)*dx1(2)*dx2(3)*dr;
    sim.disp(:,idx2(1),idx2(2),idx1(3))=sim.disp(:,idx2(1),idx2(2),idx1(3))+dx2(1)*dx2(2)*dx1(3)*dr;
    sim.disp(:,idx2(1),idx2(2),idx2(3))=sim.disp(:,idx2(1),idx2(2),idx2(3))+dx2(1)*dx2(2)*dx2(3)*dr;
end
sim.disp(:,2,:,:)=sim.disp(:,2,:,:)+sim.disp(:,n+2,:,:);
sim.disp(:,n+1,:,:)=sim.disp(:,n+1,:,:)+sim.disp(:,1,:,:);
sim.disp(:,:,2,:)=sim.disp(:,:,2,:)+sim.disp(:,:,n+2,:);
sim.disp(:,:,n+1,:)=sim.disp(:,:,n+1,:)+sim.disp(:,:,1,:);
sim.disp(:,:,:,2)=sim.disp(:,:,:,2)+sim.disp(:,:,:,n+2);
sim.disp(:,:,:,n+1)=sim.disp(:,:,:,n+1)+sim.disp(:,:,:,1);

sim.disp=sim.disp(:,2:n+1,2:n+1,2:n+1);
if 1
    figure; imagesc(reshape(mean(sim.disp(1,:,:,:),4),n,n)'); axis xy square; colorbar
    figure; imagesc(reshape(mean(sim.disp(2,:,:,:),4),n,n)'); axis xy square; colorbar
    figure; imagesc(reshape(mean(sim.disp(3,:,:,:),4),n,n)'); axis xy square; colorbar    
end
sim.div=divergence(reshape(sim.disp(1,:,:,:),[n,n,n]),...
                   reshape(sim.disp(2,:,:,:),[n,n,n]),...
                   reshape(sim.disp(3,:,:,:),[n,n,n]));
[sim.curlx,sim.curly,sim.curlz,sim.cav]=curl(reshape(sim.disp(1,:,:,:),[n,n,n]),...
                   reshape(sim.disp(2,:,:,:),[n,n,n]),...
                   reshape(sim.disp(3,:,:,:),[n,n,n]));
if 1
    figure; imagesc(-reshape(sum(sim.div(:,:,:),3),n,n)'); axis xy square; colorbar
end



end
