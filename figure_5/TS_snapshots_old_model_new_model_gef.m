%using for Figure 5A
folder = "./basegory_newreac0_gef/"
filename='dynN80gef700centered1.dat';
filenamecentroid='dynN80gef700centroid1.dat';
figurename='patch_1hr_mobility_oldmodel_gef700.pdf'
tpoints=[5 20 40 60];
nspecies=1;
save=1;


%using for Figure 5D
%folder = "./basegory_k2b_0p0625x_k5a_10x_gef_long/"
%filename='dynN80gef15centered2.dat';
%filenamecentroid='dynN80gef15centroid2.dat';
%figurename='patch_30min_mobility_1minsamp_k2b_0p0625x_k5a_10x_gef15.pdf'
%tpoints=[3 10 20 30];
%nspecies=1;
%save=1;

%using for Figure 5F 
%folder = "./basegory_k7_50_gef_long/"
%filename='dynN80gef25centered2.dat';
%filenamecentroid='dynN80gef25centroid2.dat';
%figurename='patch_30min_mobility_1minsamp_k7_50_gef25.pdf'
%tpoints=[3 10 20 30];
%nspecies=1;
%save=1;

%using for Figure 5H
%folder = "./basegory_newreac_gef/"
%filename='dynN80gef15centered1.dat';
%filenamecentroid='dynN80gef15centroid1.dat';
%figurename='patch_30min_mobility_1minsamp_newreac_gef15.pdf'
%tpoints=[3 10 20 30];
%nspecies=1;
%save=1;


params=["10000"]; 
figure('DefaultAxesFontSize',8)
nparams=size(params,2);

%full figure position
L=200;B=10;W=400;H= W/4 ;
%[loffset+w*(1+fracsep)*(i-1) b fracw/ncol h] % left, bottom, width, height.
loffset=0.02;
fracsep=0.05;
b=0.05; 
fracw=0.9;
ncol=length(tpoints);
w=fracw/ncol;
%for a squared subplot: h*H=W*w
h=w*W/H;

close all;
%could center distributions with center_patch.py
file= strcat(filename);
tseries1l=load(strcat(folder,file));
filec= strcat(filenamecentroid);
centroidspy=load(strcat(folder,filec));
centroidspy=centroidspy/10;


interval = 60;
N=sqrt(size(tseries1l,2));

Nf=size(tseries1l,1)/nspecies;
L=8;
dx=L/N
dmax=6 %maximum jump to consider
%put tseriesSL into 3 dimensional matrix form

tseries =[];
for species = 1:nspecies
    for k=1:Nf
        tseries(species,k,:,:)=reshape(tseries1l( species + (k-1)*nspecies ,:), N, N);
    end

end

%Get positions of particles from grid based simulations
maxpart = 5000;
timeanalysis=[1:Nf];
xst=-ones(maxpart,length(timeanalysis));
yst=-ones(maxpart,length(timeanalysis));
for k=1:length(timeanalysis)
    xs=[];
    ys=[];
    for i=1:N
        for j=1:N
            if tseries(1,k,i,j)>0
                xcell=(i-1)*dx + rand(1,tseries(1,k,i,j))*dx;
                ycell=(j-1)*dx + rand(1,tseries(1,k,i,j))*dx;	
                xs=[xs xcell];
                ys=[ys ycell];
            end
        end
    end
    display(k)
    %display(xs)
    xst(1:length(xs),k)=xs;
    yst(1:length(ys),k)=ys;
end	
xst(xst==-1)=nan;
yst(yst==-1)=nan;


%Create figure
zmax=5000;

%TIMESERIES OF PROTEIN DISTRIBUTION

[X,Y] = meshgrid(1:N);
X=X*dx - dx;
Y=Y*dx - dx;
centroids=[]
h_fig=figure(1);
set(h_fig,'position',[L,B,W,H]);

j=1;
for i = 1:length(tpoints)
    subplot(nparams,length(tpoints),i+(j-1)*length(tpoints))          
    %dots plot
    scatter(yst(:,tpoints(i)),xst(:,tpoints(i)),2,'red','filled')
    box on;
    set(gca,'xtick',[],'ytick',[],'LineWidth',1)
    
    xlim([0 L])
    ylim([0 L])

    tot42T = squeeze(tseries(1,tpoints(i),:,:));	
    %heat plot
    %surf(X,Y,tot42T/dx^2,'FaceColor','interp','EdgeColor', 'none');
    %xlim([0 L-dx])
    %ylim([0 L-dx])
    %zlim([0 zmax]);
    %view(0,90);
    %set(gca,'XColor', 'none','YColor','none')
    %caxis([0 2500])

    hold on
    
    %{
    %get centroid       
    tot42Tnorm=tot42T/sum(tot42T(:));
    [m,n]=size(tot42Tnorm);
    [I,J]=ndgrid(1:m,1:n);
    centroid=[dot(I(:),tot42Tnorm(:))*dx-dx,  dot(J(:),tot42Tnorm(:))*dx-dx];
    centroids=[centroids;centroid];
    
    %plot centroid
    plot3(centroids(:,2),centroids(:,1),1500*ones(size(centroids,1)),'-o','Color','black','LineWidth',1.0,'MarkerSize',2,'MarkerFaceColor','black')
    %plot3(centroids([1 end],2),centroids([1 end],1),1500*ones(2),'-o','Color','black','LineWidth',1.0,'MarkerSize',2,'MarkerFaceColor','black')
    %plot the first centroid
    plot3(centroids(1,2),centroids(1,1),1500,'o','Color','green','LineWidth',1.0,'MarkerSize',2,'MarkerFaceColor','green')
    %}
    
    %get all centroids up to this point
    if i == 1
        tot42T = squeeze(tseries(1,tpoints(i),:,:));	
        tot42Tnorm=tot42T/sum(tot42T(:));
        [m,n]=size(tot42Tnorm);
        [I,J]=ndgrid(1:m,1:n);
        %centroid1=[dot(I(:),tot42Tnorm(:))*dx-dx,  dot(J(:),tot42Tnorm(:))*dx-dx];
        %plot3(centroid1(2),centroid1(1),1500*ones(size(centroid1,1)),'-o','Color','green','LineWidth',1.0,'MarkerSize',2,'MarkerFaceColor','green')
        
    end 
    if i > 1
        for k = tpoints(i-1):tpoints(i)
            tot42T = squeeze(tseries(1,k,:,:));
            tot42Tnorm=tot42T/sum(tot42T(:));
            [m,n]=size(tot42Tnorm);
            [I,J]=ndgrid(1:m,1:n);
            centroid=[dot(I(:),tot42Tnorm(:))*dx-dx,  dot(J(:),tot42Tnorm(:))*dx-dx];
            centroids=[centroids;centroid];
            
        end
        
        %plot3(centroids(:,2),centroids(:,1),1500*ones(size(centroids,1)),'-o','Color','black','LineWidth',0.75,'MarkerSize',0.1,'MarkerFaceColor','black')
        %plot3(centroids(end,2),centroids(end,1),1500,'-o','Color','black','LineWidth',0.75,'MarkerSize',2,'MarkerFaceColor','black')
        %plot3(centroid1(2),centroid1(1),1500*ones(size(centroid1,1)),'-o','Color','green','LineWidth',0.75,'MarkerSize',2,'MarkerFaceColor','green')
        
    end
    
    if i > 1
        for k = tpoints(1):tpoints(i)-1
            if norm(centroidspy(k+1,:) - centroidspy(k,:)) < dmax
                plot3(centroidspy(k:k+1,1),centroidspy(k:k+1,2),1500*ones(2),'-o','Color','black','LineWidth',0.75,'MarkerSize',0.1,'MarkerFaceColor','black')
            end
        end
    end
    %plot3(centroidspy(tpoints(1):tpoints(i),1),centroidspy(tpoints(1):tpoints(i),2),1500*ones(size(tpoints(1):tpoints(i))),'-o','Color','black','LineWidth',0.75,'MarkerSize',0.1,'MarkerFaceColor','black')
    
    plot3(centroidspy(tpoints(i),1),centroidspy(tpoints(i),2),1500,'-o','Color','black','LineWidth',0.75,'MarkerSize',2,'MarkerFaceColor','black')
    plot3(centroidspy(tpoints(1),1),centroidspy(tpoints(1),2),1500,'-o','Color','green','LineWidth',0.75,'MarkerSize',2,'MarkerFaceColor','green')

    
    hold off
    %colorbar
    %title({[num2str(tpoints(i)) ' min']},'FontSize',8)
    %title({[num2str(sum(tot42T(:))) ' Tot Cdc42T']},'FontSize',6)
    h_subplot(i+(j-1)*length(tpoints)) = subplot(nparams,length(tpoints),i+(j-1)*length(tpoints));
    set(h_subplot(i+(j-1)*length(tpoints)), 'position', [loffset+w*(1+fracsep)*(i-1) b w h] );
    [loffset+w*(1+fracsep)*(i-1) b w*W h*H];
    %get(h_subplot(i+(j-1)*length(tpoints)), 'position' )

end

if save ==1 
    saveas(h_fig,figurename)
end
