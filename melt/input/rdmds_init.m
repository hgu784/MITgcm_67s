%This program creates open boundary prescription files for the PIG
%experiment based on the 10 yrs spin-up run using OBCS with
%U,V = 0, Tref, Sref

%following variables
%-----3D fields-----
% T Temperature (C)
% S Salinity (psu)
% U u-velocity (m/s)
% PH ocean pressure (or atm geopotential)


addpath /home/toshiki/MITgcm/utils/matlab/cs_grid/read_cs/
addpath /home/toshiki/MITgcm/utils/matlab/cs_grid/
addpath /home/toshiki/MITgcm/utils/matlab/
%% seasonal
clear
nx = 60;    delx = 1;   X = nx*delx;
ny = 100;    dely = 1;   Y = ny*dely;
nz = 55;    delz = 20;  Z = nz*delz;

tt = 12*3*60;

%x axis
x = 1:delx:X;
%y axis
y = 1:dely:Y;
%z axis [m]
z = delz:delz:Z;


T_OBC(1:nx,1:nz) = NaN;


T_sfc = -1;
T_bottom = 1.2;
del_T = (T_bottom - T_sfc)/(400);

for iz = 1:nz
    tref(iz) = T_sfc + del_T*((iz-15)*delz);
    if (iz)<=15
        tref(iz)=-1;
    end
    if iz>=35
        tref(iz) =1.2;
    end
end



for iz = 1:nz
    T_OBC(1:nx,iz) = tref(iz);
end


%  Linear S profile
%  From Smin at sfc increasing to Smax at bottom
%  Take Smin and Smax from the OBC of the original PIG experiment

% Initialize variable
S_OBC(1:nx,1:nz) = NaN;

S_sfc = 34;
S_bottom = 34.7;
del_S = (S_bottom - S_sfc)/(400);

for iz = 1:nz
    sref(iz) = S_sfc + del_S*((iz-15)*delz);
    if iz<=15
        sref(iz)=34;
    end
    if iz>=35
        sref(iz) =34.7;
    end
end

for iz = 1:nz
    S_OBC(1:nx,iz) = sref(iz);
end

To = zeros(nx,nz, tt);
So = zeros(nx,nz, tt);


for it = 1:tt
    To(:,:,it) = T_OBC;
    So(:,:,it) = S_OBC;
end

fid=fopen('theta_2.obw','w','b');  fwrite(fid,To,'real*8'); fclose(fid);
fid=fopen('salt_2.obw','w','b');  fwrite(fid,So,'real*8'); fclose(fid);




%%



%Set the grid size;
nx = 60;    delx = 1;   X = nx*delx;
ny = 100;    dely = 1;   Y = ny*dely;
nz = 55;    delz = 20;  Z = nz*delz;

%x axis
x = 1:delx:X;
%y axis
y = 1:dely:Y;
%z axis [m]
z = delz:delz:Z;


%Note:
%Ensure that the volume flux is zero at the open boundary


%% ________________________________________________________________________
%  First try: linear velocity profile.
%  From U = -0.025 ms-1 (outwards) to  0.025 ms-1
%  Constant in y direction (meridionally)
%  Do not take exactly Umin and Umax from the OBC of the original
%  experiment: Take modified numbers to that Umin = -Umax and volume is conserved

% Initialize variable
V_OBC(1:nx,1:nz) = NaN;


v_sfc = 0.025;                         %U at sfc
v_bottom = -0.025;                       %U at bottom
del_v = (v_bottom - v_sfc)/((nz-1)*delz);   %delu/delz

for iz = 1:nz
    V_OBC(1:nx,iz) = v_sfc + del_v*((iz-1)*delz);
end


%  Linear T profile
%  From Tmin at sfc increasing to Tmax at bottom
%  Take Tmin and Tmax from the OBC of the original PIG experiment

% Initialize variable
T_OBC(1:nx,1:nz) = NaN;


T_sfc = -1;
T_bottom = 1.2;
del_T = (T_bottom - T_sfc)/(400);

for iz = 1:nz
    
    
    tref(iz) = T_sfc + del_T*((iz-15)*delz);
    if (iz)<=15
        tref(iz)=-1;
    end
    if iz>=35
        tref(iz) =1.2;
    end
end



for iz = 1:nz
    T_OBC(1:nx,iz) = tref(iz);
end


%  Linear S profile
%  From Smin at sfc increasing to Smax at bottom
%  Take Smin and Smax from the OBC of the original PIG experiment

% Initialize variable
S_OBC(1:nx,1:nz) = NaN;

S_sfc = 34;
S_bottom = 34.7;
del_S = (S_bottom - S_sfc)/(400);

for iz = 1:nz
    
    
    sref(iz) = S_sfc + del_S*((iz-15)*delz);
    if iz<=15
        sref(iz)=34;
    end
    if iz>=35
        sref(iz) =34.7;
    end
end

for iz = 1:nz
    S_OBC(1:nx,iz) = sref(iz);
end




%% Initial T & S conditions

% Take western open boundary conditions for T & S
% and assume no change in x-direction

T_init = zeros(nx,ny,nz);
S_init = zeros(nx,ny,nz);

for iy = 1:ny
    T_init(:,iy,:) = T_OBC;
    S_init(:,iy,:) = S_OBC;
end



%% Print OBCS files for T, S

% fid=fopen('uvel.obw','w','b');  fprintf(fid,'%10.4f',U_OBC);fclose(fid);
fid=fopen('vvel.obw','w','b');  fwrite(fid,V_OBC,'real*8'); fclose(fid);
% fid=fopen('uvel.obw','w','b');  fwrite(fid,U_OBC,'real*8'); fclose(fid);
fid=fopen('theta.obw','w','b');  fwrite(fid,T_OBC,'real*8'); fclose(fid);
fid=fopen('salt.obw','w','b');  fwrite(fid,S_OBC,'real*8'); fclose(fid);
% fid=fopen('theta.obw','w');  fprintf(fid,'%10.4f',T_OBC);fclose(fid);
% fid=fopen('salt.obw','w');  fprintf(fid,'%10.4f',S_OBC);fclose(fid);



%% Print init files for T, S

fid=fopen('theta.init','w','b');fwrite(fid,T_init,'real*8');fclose(fid);
fid=fopen('salt.init','w','b');fwrite(fid,S_init,'real*8');fclose(fid);


% for iz = 1:nz;
%     fprintf(fid_T,'%10.4f',T_init(:,:,iz));
%     fprintf(fid_S,'%10.4f',S_init(:,:,iz));
% end
%
% fclose(fid_T);
% fclose(fid_S);
%%
a = readbin('/home/toshiki/MITgcm/run/jim_run3/input/theta.obw', [60 55], 1, 'real*8');
%%
figure
quiver()