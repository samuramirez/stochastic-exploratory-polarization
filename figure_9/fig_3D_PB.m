
filename='fig_3D_transition/3d_newreac_gef100_k4a_0p01_04.xyz';
figurename='3d_newreac_gef100_k4a_0p01.pdf'

filename='fig_3D_transition/3d_newreac_gef100_k4a_0p002_04.xyz';
figurename='3d_newreac_gef100_k4a_0p002.pdf'
filename='fig_3D_transition/3d_newreac_gef100_k4a_0p001_04.xyz';
figurename='3d_newreac_gef100_k4a_0p001.pdf'

filename='fig_3D_transition/3d_newreac_gef100_k4a_0p02_01.xyz';
figurename='3d_newreac_gef100_k4a_0p02.pdf'

filename='fig_3D_transition/3d_newreac_gef100_k4a_0p005_04.xyz';
figurename='3d_newreac_gef100_k4a_0p005.pdf'



nframes=60; %not the actual frames, but used to define containers in read_molPos3

[t,positions]=read_molPos3(filename,nframes);

L=200;B=200;W=600;H= W/4 ;

hf=figure(1);
set(hf,'position',[L,B,W,H]);

snapshots=[3,10,17,23];

for j=1:length(snapshots)
    i=snapshots(j);
    subplot(1,4,j)
    
    % Code to find active Cdc42, may need to modify
    if ~isempty(positions.Cdc42T{i}) && ~isempty(positions.BemGEF42{i})
        all_active_x = [positions.Cdc42T{i}(:,1);positions.BemGEF42{i}(:,1)];
        all_active_y = [positions.Cdc42T{i}(:,2);positions.BemGEF42{i}(:,2)];
        all_active_z = [positions.Cdc42T{i}(:,3);positions.BemGEF42{i}(:,3)];
    elseif isempty(positions.Cdc42T{i}) && ~isempty(positions.BemGEF42{i})
        all_active_x = positions.BemGEF42{i}(:,1);
        all_active_y = positions.BemGEF42{i}(:,2);
        all_active_z = positions.BemGEF42{i}(:,3);
    elseif ~isempty(positions.Cdc42T{i}) && isempty(positions.BemGEF42{i})
        all_active_x = positions.Cdc42T{i}(:,1);
        all_active_y = positions.Cdc42T{i}(:,2);
        all_active_z = positions.Cdc42T{i}(:,3);
    else
        all_active_x = nan;
        all_active_y = nan;
        all_active_z = nan;
        
    end
%     vis3dhist(all_active_x,all_active_y,all_active_z,2.5,2.5,2.5,'active Cdc42 molecules');
%     axis image;
%     axis([0 3.6340 0 3.6340]); hold on;
%     pbaspect([1 .5 .8518])
%     set(gca,'zlim',[0 40]);
    r = 2.2567;
    boxLim= 2.5;
    X=all_active_x;
    Y=all_active_y;
    Z=all_active_z;
    X=X-r; Y=Y-r; Z=Z-r;
    [x,y,z]=sphere(300);
    x=x*r; y=y*r; z=z*r;

    surf(x,y,z,'facecolor',[1 1 1],'FaceAlpha', 0.3,'FaceLighting','gouraud','edgecolor','none'); 
    camlight 
    hold on;
    plot3(X,Y,Z,'.','Color','red')
    grid off; box off; axis equal;
    view(0,0);
    hold off;
    %title('view 1')
    set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[]);%,'clim',[0 40]);
    %axis(boxLim*[-1 1 -1 1 -1 1]);
    axis off;
    %text(-15,-3.5,-5,sprintf('%i seconds',t(i)),'fontsize',14);
%     print(gcf,'-dpng',sprintf('%s/frame%03i.png',moviedir,i),'-r100');
%     pause(0.1);

    %drawnow;
    %j=j+1;
    

end
%saveas(hf,figurename)
