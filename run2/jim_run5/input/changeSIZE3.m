%% [changesize3] -> [rdmds_init] -> [gendata]
clear
clc

nx = 120;
ny = 200;

% H = readbin('/home/toshiki/MITgcm/run/jim_run_original/input/Hinit1100.box', [60 160], 1, 'real*8');
H = rdmds('/home/toshiki/MITgcm/run2/jim_run5/input/surfDiag.0018921600');
H = H(:,:,6);
h = H(2:59, 2:60);
J = imresize(h,[nx-2 nx-1]); % (60*60 to 120*120)              　% 注意↓
one = ones(nx-2, 1)*400; one2 = ones(1, nx)*400; zero = zeros(nx, ny-nx);
H2 = horzcat(one, J); H2 = vertcat(one2, H2, one2); H2 = horzcat(H2, zero);
fid = fopen('Hinit60y.box','w','b'); fwrite(fid,H2,'real*8'); fclose(fid);

                                  % 注意↓
% hm1 = ones(nx, nx); hm0 = zeros(nx, ny-nx);
% HM2 = horzcat(hm1, hm0); HM2(:, 1) = -1; HM2([1 end], :) = -1;
% fid = fopen('HMASK60.box','w','b'); fwrite(fid,HM2,'real*8'); fclose(fid);

                                  % 注意↓
im1 = ones(nx, nx); im0 = zeros(nx, ny-nx);
IM2 = horzcat(im1, im0);
fid = fopen('IM_CALVE60.box','w','b'); fwrite(fid,IM2,'real*8'); fclose(fid);

clearvars -except nx ny H2 HM2 IM2 