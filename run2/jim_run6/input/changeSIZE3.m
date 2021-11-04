%% [changesize3] -> [rdmds_init] -> [gendata]
clear
clc

nx = 240;
ny = 400;

% H = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/Hinit1100.box', [60 160], 1, 'real*8');
H = rdmds('/home/toshiki/MITgcm/run2/jim_run6/input/surfDiag.0018921600');
H = H(:,:,6);
h = H(2:59, 2:60);
J = imresize(h,[nx-2 nx-1]);                                   % 注意↓
one = ones(nx-2, 1)*400; one2 = ones(1, nx)*400; zero = zeros(nx, ny-nx);
H2 = horzcat(one, J); H2 = vertcat(one2, H2, one2); H2 = horzcat(H2, zero);
fid = fopen('Hinit60y.box','w','b'); fwrite(fid,H2,'real*8'); fclose(fid);

%                                   % 注意↓
% hm1 = ones(nx, nx); hm0 = zeros(nx, ny-nx);
% HM2 = horzcat(hm1, hm0); HM2(:, 1) = -1; HM2([1 end], :) = -1;
% fid = fopen('HMASK60.box','w','b'); fwrite(fid,HM2,'real*8'); fclose(fid);

                                  % 注意↓
im1 = ones(nx, nx); im0 = zeros(nx, ny-nx);
IM2 = horzcat(im1, im0);
fid = fopen('IM_CALVE60.box','w','b'); fwrite(fid,IM2,'real*8'); fclose(fid);

clearvars -except nx ny H2 HM2 IM2 

%%
clear
% 460 600 130

a = readbin('/home/toshiki/MITgcm/run2/pig/latlon_run23/Bedrock_BedMachine_ext_qsg_yosi3_init.bin',[460 600],1,'real*4');
b = imresize(a,[115 150]);

mp(a)
figure
mp(b)




