a = readbin('initT_latlon_qsg.bin', [460 600 130],1,'real*4');

%%
% clear
%
% t = readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/initT_latlon_qsg.bin',[460 600 130],1,'real*4');
% fid = fopen('theta.init','w','b'); fwrite(fid,t,'real*4'); fclose(fid);
% s = readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/initS_latlon_qsg.bin',[460 600 130],1,'real*4');
% fid = fopen('salt.init','w','b'); fwrite(fid,s,'real*4'); fclose(fid);


%%
clear
t = readbin('Z3D_S2.bin',[460 130 1440],1,'real*4');

find(t)




%%

clear

OBNsFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBNs_200_for_cilan_rev.bin',[460 130 1440],1,'real*4');Ns = mean(OBNsFile,3);
OBNtFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBNt_200_for_cilan_rev.bin',[460 130 1440],1,'real*4');Nt = mean(OBNtFile,3);
OBNuFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBNu_200_for_cilan_rev.bin',[460 130 1440],1,'real*4');Nu = mean(OBNuFile,3);
OBNvFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBNv_200_for_cilan_rev.bin',[460 130 1440],1,'real*4');Nv = mean(OBNvFile,3);

OBWsFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBEs_200_for_cilan_rev.bin',[600 130 1440],1,'real*4');Ws = mean(OBWsFile,3);
OBWtFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBEt_200_for_cilan_rev.bin',[600 130 1440],1,'real*4');Wt = mean(OBWtFile,3);
OBWuFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBEu_200_for_cilan_rev.bin',[600 130 1440],1,'real*4');Wu = mean(OBWuFile,3);
OBWvFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBEv_200_for_cilan_rev.bin',[600 130 1440],1,'real*4');Wv = mean(OBWvFile,3);

%%
ns = imresize(Ns, [115 130],'nearest'); ws = imresize(Ws, [150 130],'nearest');
nt = imresize(Nt, [115 130],'nearest'); wt = imresize(Wt, [150 130],'nearest');
nu = imresize(Nu, [115 130],'nearest'); wu = imresize(Wu, [150 130],'nearest');
nv = imresize(Nv, [115 130],'nearest'); wv = imresize(Wv, [150 130],'nearest');

nx = 115; ny = 150; nz = 130; del_t = 1440;

ns1 = zeros(nx,nz,del_t); ws1 = zeros(ny,nz,del_t);
nt1 = zeros(nx,nz,del_t); wt1 = zeros(ny,nz,del_t);
nu1 = zeros(nx,nz,del_t); wu1 = zeros(ny,nz,del_t);
nv1 = zeros(nx,nz,del_t); wv1 = zeros(ny,nz,del_t);


for i = 1:del_t
    ns1(:,:,i) = ns; ws1(:,:,i) = ws;
    nt1(:,:,i) = nt; wt1(:,:,i) = wt;
    nu1(:,:,i) = nu; wu1(:,:,i) = wu;
    nv1(:,:,i) = nv; wv1(:,:,i) = wv;
end

%%
Z115 = zeros(115,130,1440); Z150 = zeros(150,130, 1440);
%%
fid=fopen('FOBNs_200_for_cilan_rev2.bin','w','b');fwrite(fid,ns1,'real*4');fclose(fid);
fid=fopen('FOBNt_200_for_cilan_rev2.bin','w','b');fwrite(fid,nt1,'real*4');fclose(fid);
fid=fopen('FOBNu_200_for_cilan_rev2.bin','w','b');fwrite(fid,nu1,'real*4');fclose(fid);
fid=fopen('FOBNv_200_for_cilan_rev2.bin','w','b');fwrite(fid,nv1,'real*4');fclose(fid);

fid=fopen('FOBWs_200_for_cilan_rev2.bin','w','b');fwrite(fid,ws1,'real*4');fclose(fid);
fid=fopen('FOBWt_200_for_cilan_rev2.bin','w','b');fwrite(fid,wt1,'real*4');fclose(fid);
fid=fopen('FOBWu_200_for_cilan_rev2.bin','w','b');fwrite(fid,wu1,'real*4');fclose(fid);
fid=fopen('FOBWv_200_for_cilan_rev2.bin','w','b');fwrite(fid,wv1,'real*4');fclose(fid);

fid=fopen('Z3D_SS2.bin','w','b');fwrite(fid,Z115,'real*4');fclose(fid);
fid=fopen('Z3D_EE2.bin','w','b');fwrite(fid,Z150,'real*4');fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Efile
clear
OBEsFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBEs_600x130x1440_new.bin',[600 130 1440],1,'real*4');Es = mean(OBEsFile,3);
OBEtFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBEt_600x130x1440_2_new_yosi2.bin',[600 130 1440],1,'real*4');Et = mean(OBEtFile,3);
OBEuFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBEu_600x130x1440_2_new_yosi2.bin',[600 130 1440],1,'real*4');Eu = mean(OBEuFile,3);
OBEvFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/FOBEv_600x130x1440_new.bin',[600 130 1440],1,'real*4');Ev = mean(OBEvFile,3);

es = imresize(Es, [150 130],'nearest');
et = imresize(Et, [150 130],'nearest');
eu = imresize(Eu, [150 130],'nearest');
ev = imresize(Ev, [150 130],'nearest');

nx = 115; ny = 150; nz = 130; del_t = 1440;

es1 = zeros(ny,nz,del_t);
et1 = zeros(ny,nz,del_t);
eu1 = zeros(ny,nz,del_t);
ev1 = zeros(ny,nz,del_t);

for i = 1:del_t
    es1(:,:,i) = es;
    et1(:,:,i) = et;
    eu1(:,:,i) = eu;
    ev1(:,:,i) = ev;
end

fid=fopen('FOBEs_200_for_cilan_rev2.bin','w','b');fwrite(fid,es1,'real*4');fclose(fid);
fid=fopen('FOBEt_200_for_cilan_rev2.bin','w','b');fwrite(fid,et1,'real*4');fclose(fid);
fid=fopen('FOBEu_200_for_cilan_rev2.bin','w','b');fwrite(fid,eu1,'real*4');fclose(fid);
fid=fopen('FOBEv_200_for_cilan_rev2.bin','w','b');fwrite(fid,ev1,'real*4');fclose(fid);



%%
clear

OBSsFile=readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/Z3D_S2.bin',[460 130 1440],1,'real*4');Es = mean(OBSsFile,3);







