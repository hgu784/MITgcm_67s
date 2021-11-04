%%
clear 
close
%%%%%%% 解像度あげる %%%%%%%%%
%%%%%%% ~/jim_runOO/input/ の中で %%%%%%%%%%%%%%%%%%%%%%%%
H = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/Hinit1100.box', [60 160], 1, 'real*8');
h = H(2:59, 2:60);
J = imresize(h,[118 119]); % (60*60 to 120*120)
one = ones(118, 1)*400; one2 = ones(1, 120)*400; zero = zeros(120, 80);
H2 = horzcat(one, J); H2 = vertcat(one2, H2, one2); H2 = horzcat(H2, zero);
fid = fopen('Hinit1100.box','w','b'); fwrite(fid,H2,'real*8'); fclose(fid);

% B = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/BATHY60.box', [60 160], 1, 'real*8');
B2 = ones(120, 200)*-1100;
B2(:, 1) = 0; B2(1, :) = 0; B2(end, :) = 0; 
fid = fopen('BATHY60.box','w','b'); fwrite(fid,B2,'real*8'); fclose(fid);

% UF = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/UFACEMASK60.box', [60 160], 1, 'real*8');
UF2 = ones(120, 200)*-1;
fid = fopen('UFACEMASK60.box','w','b'); fwrite(fid,UF2,'real*8'); fclose(fid);

% VF = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/VFACEMASK60.box', [60 160], 1, 'real*8');
VF2 = UF2; VF2(:, 2) = 3;
fid = fopen('VFACEMASK60.box','w','b'); fwrite(fid,VF2,'real*8'); fclose(fid);

VD = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/VDIRICH60.box', [60 160], 1, 'real*8');
vd = VD(:, 2); vd = imresize(vd, [120 1]);
VD2 = zeros(120, 200); VD2(:,2) = vd;
fid = fopen('VDIRICH60.box','w','b'); fwrite(fid,VD2,'real*8'); fclose(fid);

% HBC = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/HBCy60.box', [60 160], 1, 'real*8');
HBC2 = zeros(120, 200); HBC2(:, 2) = 1200;
fid = fopen('HBCy60.box','w','b'); fwrite(fid,HBC2,'real*8'); fclose(fid);

% HM = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/HMASK60.box', [60 160], 1, 'real*8');
hm1 = ones(120, 120); hm0 = zeros(120, 80);
HM2 = horzcat(hm1, hm0); HM2(:, 1) = -1; HM2(1, :) = -1; HM2(end, :) = -1;
fid = fopen('HMASK60.box','w','b'); fwrite(fid,HM2,'real*8'); fclose(fid);

e0 = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/etainit0.round.bin', [60 160], 1, 'real*8');
e = e0(2:59, 2:60);
K = imresize(e,[118 119]); % (60*60 to 120*120)
one = ones(118, 1)*3.2688; one2 = ones(1, 120)*3.2688; zero = zeros(120, 80);
e02 = horzcat(one, K); e02 = vertcat(one2, e02, one2); e02 = horzcat(e02, zero);
fid = fopen('etainit0.round.bin','w','b'); fwrite(fid,e02,'real*8'); fclose(fid);

topo = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/shelftopo.round.bin', [60 160], 1, 'real*8');
tp = topo(2:59, 2:60);
L = imresize(tp,[118 119], 'box');
one = ones(118, 1)*-360; one2 = ones(1, 120)*-360; zero = zeros(120, 80);
topo2 = horzcat(one, L); topo2 = vertcat(one2, topo2, one2); topo2 = horzcat(topo2, zero);
fid = fopen('shelftopo.round.bin','w','b'); fwrite(fid,topo2,'real*8'); fclose(fid);

eta = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/etainit.round.bin', [60 160], 1, 'real*8');
et = eta(2:59, 2:60);
K = imresize(et,[118 119],'bilinear'); % (60*60 to 120*120)
one = ones(118, 1)*3.2688; one2 = ones(1, 120)*3.2688; zero = zeros(120, 80);
eta2 = horzcat(one, K); eta2 = vertcat(one2, eta2, one2); eta2 = horzcat(eta2, zero);
fid = fopen('etainit.round.bin','w','b'); fwrite(fid,eta2,'real*8'); fclose(fid);

mass = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/shelficemassinit.bin', [60 160], 1, 'real*8');
ma = mass(2:59, 2:60);
I = imresize(ma,[118 119]); % (60*60 to 120*120)
one = ones(118, 1)*366800; one2 = ones(1, 120)*366800; zero = zeros(120, 80);
mass2 = horzcat(one, I); mass2 = vertcat(one2, mass2, one2); mass2 = horzcat(mass2, zero);
fid = fopen('shelficemassinit.bin','w','b'); fwrite(fid,mass2,'real*8'); fclose(fid);

% IM = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/IM_CALVE60.box', [60 160], 1, 'real*8');
im1 = ones(120, 120); im0 = zeros(120, 80);
IM2 = horzcat(im1, im0);
fid = fopen('IM_CALVE60.box','w','b'); fwrite(fid,IM2,'real*8'); fclose(fid);

clearvars -except H2 B2 UF2 VF2 VD2 HBC2 HM2 e02 topo2 eta2 mass2 IM2








