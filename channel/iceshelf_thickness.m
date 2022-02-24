%% video

clear
folder = '~/PATH/';    % <-- change
file = dir(fullfile(folder,'surfDiag.*.data'));
str = strcat(folder,'surfDiag');

for i = 1 :36: max(size(file))
    
    f = figure(1);
    f.WindowState = 'maximized';
    
    a = nan(60,60);
    surf = rdmds(str, str2double(file(i).name(11:end-5))); 
    s = surf(2:59,2:60,6); a(2:59,2:60) = s;

    ss = mesh(-a,'FaceAlpha','0.5'); ss.FaceColor = 'flat';colorbar
    zlim([-1200 0]); caxis([-900 -100]); colormap(jet) ;
    ax = gca; ax.YDir = 'reverse'; ax.FontSize = 18; ax.TickLength = [0 0];
    xticks([0 20 40 60]); yticks([0 20 40 60]); 
    xlabel('Along shelf (km)');ylabel('Across shelf (km)');zlabel('Depth (m)');
    
    title([num2str(i*10), 'days'])
    
end