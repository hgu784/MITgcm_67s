
%Verion of gendata.m modified by DNG
%This is a matlab script that generates the input data


% the configuation approximately the ISOMIP experiment no. 1
% require matlab functions for equation of state
clear

% Dimensions of grid
nx=60; 
ny=100;
nz=55;
delz = 20;

hfacMin = 0.2;
%mwct = 3;

dlat = 0.125/64*4; dy=dlat;
dlon = 0.015625; dx=dlon;

%eos = 'linear';
eos = 'jmd95z';
% eos = 'mdjwf';

acc = 'real*8';

long = [-105.5:dlon:-105.5];
lonc = long+dlon/2;
latg = [-75.4457:dlat:-73.8809-dlat];
latc = latg+dlat/2;
size(latc);



dz = delz*ones(1,nz);
zgp1 = [0,cumsum(dz)];
zc = .5*(zgp1(1:end-1)+zgp1(2:end));
zg = zgp1(1:end-1);
dz = diff(zgp1);
% sprintf('delZ = %d * %7.6g,',nz,dz)

%%%%%%%%% STREAMICE FILES %%%%%%%%%%%%%%%

fid = fopen('Hinit1100.box','r','b'); H_streamice = fread(fid,[nx ny],'real*8'); fclose(fid);

bathy = -1100 * ones(ny,nx); 
bathy(1,:) = 0; bathy(:,[1 end]) = 0;
fid = fopen('BATHY60.box','w','b'); fwrite(fid,bathy','real*8'); fclose(fid);

ufacemask = -1 * ones(ny,nx);
vfacemask = -1 * ones(ny,nx); vfacemask(2,:) = 3;
fid = fopen('UFACEMASK60.box','w','b'); fwrite(fid,ufacemask','real*8'); fclose(fid);
fid = fopen('VFACEMASK60.box','w','b'); fwrite(fid,vfacemask','real*8'); fclose(fid);

avgV = 1000;
x = dx/2:dx:(nx*dx-dx/2);
xmid = nx*dx/2;
Vin = (2*avgV/xmid^2) * (xmid^2-(x - xmid).^2);
vdirich = zeros(ny,nx);
vdirich(2,:) = Vin;
fid = fopen('VDIRICH60.box','w','b'); fwrite(fid,vdirich','real*8'); fclose(fid);

Hin = 1200;
Hbc = zeros(ny,nx);
Hbc(2,:) = Hin;
fid = fopen('HBCy60.box','w','b'); fwrite(fid,Hbc','real*8'); fclose(fid);

hmask = ones(ny,nx);
hmask(60+1:end,:) = 0;
hmask(1,:) = -1;
hmask(:,[1 end]) = -1;
fid = fopen('HMASK60.box','w','b'); fwrite(fid,hmask','real*8'); fclose(fid);

%%%%%%%%% stratification %%%%%%%%%%%%%%%%

fid = fopen('theta.init','r','b'); q = fread(fid,inf,'real*8'); fclose(fid); Tinit=reshape(q,[nx ny nz]);
fid = fopen('salt.init','r','b'); q = fread(fid,inf,'real*8'); fclose(fid); Sinit=reshape(q,[nx ny nz]);
sref = squeeze(Sinit(1,1,:));
tref = squeeze(Tinit(1,1,:));

%%%%%%%%%%% density %%%%%%%%%%%%%%%%%

% Gravity
gravity=9.81;
rhoConst = 1000;
k=1;
dzm = abs([zg(1)-zc(1) .5*diff(zc)]);
dzp = abs([.5*diff(zc) zc(end)-zg(end)]);
p = abs(zc)*gravity*rhoConst*1e-4;
dp = p;
kp = 0;

Rho = zeros(nz,1);

while rms(dp) > 1e-13
  phiHydF(k) = 0;
  p0 = p;
  kp = kp+1;
  for k = 1:nz
    switch eos
     case 'linear'

     case 'jmd95z'
      drho = densjmd95(sref(k),tref(k),p(k))-rhoConst;
     case 'mdjwf'
      drho = densmdjwf(sref(k),tref(k),p(k))-rhoConst;
     otherwise
      error(sprintf('unknown EOS: %s',eos))
    end
    Rho(k) = drho+rhoConst;
    phiHydC(k)   = phiHydF(k) + dzm(k)*gravity*drho/rhoConst;
    phiHydF(k+1) = phiHydC(k) + dzp(k)*gravity*drho/rhoConst;
  end
  switch eos
   case 'mdjwf'
    p = (gravity*rhoConst*abs(zc) + phiHydC*rhoConst)/gravity/rhoConst;
  end
  dp = p-p0;
end

%shelficemass = binread('HINIT60.bin',nx,ny) * 917;


fid = fopen('Hinit1100.box','r','b');
VV = fread(fid,inf,'real*8');
fclose(fid);

VV=reshape(VV, [nx ny]);

shelficemass=VV*917;
%phi0surf = zeros(nx,ny);
topo = zeros(nx,ny);

for ix=1:nx
  for iy=1:ny
%    k=max(find(abs(zg)<abs(topo(ix,iy))));
%    if isempty(k)
%      k=0;
%    end
%    if k>0
      
%      dr = -zg(k) - topo(ix,iy);
      
%      if (dr<delz/2)
%          phi0surf(ix,iy) = phiHydF(k) + (delz/2-dr) * (phiHydC(k)-phiHydF(k))/(delz/2);
%      else
%          phi0surf(ix,iy) = phiHydC(k) + (dr-delz/2) * (phiHydF(k+1)-phiHydC(k))/(delz/2);
%      end
      
%    end

     mass = shelficemass (ix,iy);
     massFuncC = rhoConst * (phiHydC/gravity + zc);
     massFuncF = rhoConst * (phiHydF/gravity + zgp1);

     k = max (find ( massFuncF < mass ));
     if (isempty(k))
         k=0;
     end
     if (k>0)
     if (mass < massFuncC(k))
      topo(ix,iy) = -zg(k) - (mass-massFuncF(k)) * delz/2 / (massFuncC(k)-massFuncF(k));
     else
      topo(ix,iy) = -zc(k) - (mass-massFuncC(k)) * delz/2 / (massFuncF(k+1)-massFuncC(k));
     end
     end

  end
end

%mass = rhoConst * (phi0surf / gravity - topo);
%mass(:,1:100) = mass(:,1:100) + 917 * 10;
%topo(:,1:100) = bathy(:,1:100)+mwct;

etainit = zeros(size(topo));

% new topography: icetopo rounded to the nearest k * deltaZ
%                 eta_init set to make difference

icetopo2 = topo;

for ix=1:nx
  for iy=1:ny
    k=max(find(abs(zg)<abs(icetopo2(ix,iy))));
    if isempty(k)
      k=0;
    else
      
      dr = 1-(-zg(k) - icetopo2(ix,iy))/delz;
      if (dr > .25)
          % bring Ro_surf *up* to closest grid face & make etainit negative
          % to compensate
          icetopo2(ix,iy) = -zg(k);
          etainit(ix,iy) = (dr-1)*delz;
      else
          % bring Ro_surf *down* to closest grid face & make etainit pos
          % to compensate
          icetopo2(ix,iy) = -zg(k+1);
          etainit(ix,iy) = (dr)*delz;
      end
       
    end
  end
end

%etainit(:,:)=0;
%icetopo2(:,1)=0;
fid = fopen('etainit0.round.bin','w','b'); fwrite(fid,etainit,'real*8'); fclose(fid);
fid = fopen('shelftopo.round.bin','w','b'); fwrite(fid,icetopo2,'real*8'); fclose(fid);
fid = fopen('etainit.round.bin','w','b'); fwrite(fid,etainit,'real*8'); fclose(fid);
fid = fopen('shelficemassinit.bin','w','b'); fwrite(fid,shelficemass,'real*8'); fclose(fid);
%fid = fopen('bathy_step.bin','w','b'); fwrite(fid,bathy,'real*8'); fclose(fid);

