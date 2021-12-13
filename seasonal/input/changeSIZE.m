clear
H = rdmds('/home/toshiki/MITgcm/run/jim_run3/run/surfDiag.0018921600');
H2 = H(:,:,6);
fid = fopen('Hinit60y.box','w','b'); fwrite(fid,H2,'real*8'); fclose(fid);






%% 
%%%%%%% ~/jim_runOO/input/ の中で %%%%%%%%%%%%%%%%%%%%%%%%
H = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/Hinit1100.box', [60 160], 1, 'real*8');
H(:, end-59:end) = [];
fid = fopen('Hinit1100.box','w','b'); fwrite(fid,H,'real*8'); fclose(fid);

B = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/BATHY60.box', [60 160], 1, 'real*8');
B(:, end-59:end) = [];
fid = fopen('BATHY60.box','w','b'); fwrite(fid,B,'real*8'); fclose(fid);

UF = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/UFACEMASK60.box', [60 160], 1, 'real*8');
UF(:, end-59:end) = [];
fid = fopen('UFACEMASK60.box','w','b'); fwrite(fid,UF,'real*8'); fclose(fid);

VF = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/VFACEMASK60.box', [60 160], 1, 'real*8');
VF(:, end-59:end) = [];
fid = fopen('VFACEMASK60.box','w','b'); fwrite(fid,VF,'real*8'); fclose(fid);

VD = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/VDIRICH60.box', [60 160], 1, 'real*8');
VD(:, end-59:end) = [];
fid = fopen('VDIRICH60.box','w','b'); fwrite(fid,VD,'real*8'); fclose(fid);

HBC = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/HBCy60.box', [60 160], 1, 'real*8');
HBC(:, end-59:end) = [];
fid = fopen('HBCy60.box','w','b'); fwrite(fid,HBC,'real*8'); fclose(fid);

HM = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/HMASK60.box', [60 160], 1, 'real*8');
HM(:, end-59:end) = [];
fid = fopen('HMASK60.box','w','b'); fwrite(fid,HM,'real*8'); fclose(fid);

e0 = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/etainit0.round.bin', [60 160], 1, 'real*8');
e0(:, end-59:end) = [];
fid = fopen('etainit0.round.bin','w','b'); fwrite(fid,e0,'real*8'); fclose(fid);

topo = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/shelftopo.round.bin', [60 160], 1, 'real*8');
topo(:, end-59:end) = [];
fid = fopen('shelftopo.round.bin','w','b'); fwrite(fid,topo,'real*8'); fclose(fid);

eta = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/etainit.round.bin', [60 160], 1, 'real*8');
eta(:, end-59:end) = [];
fid = fopen('etainit.round.bin','w','b'); fwrite(fid,eta,'real*8'); fclose(fid);

mass = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/shelficemassinit.bin', [60 160], 1, 'real*8');
mass(:, end-59:end) = [];
fid = fopen('shelficemassinit.bin','w','b'); fwrite(fid,mass,'real*8'); fclose(fid);

IM = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/IM_CALVE60.box', [60 160], 1, 'real*8');
IM(:, end-59:end) = [];
fid = fopen('IM_CALVE60.box','w','b'); fwrite(fid,IM,'real*8'); fclose(fid);


