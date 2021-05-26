
folder="./basegorymeth_rho20nmgef700/"

%THE FOLLOWING TWO FIGURES ARE PANELS D AND E OF FIGURE 2  

file= 'basegorykhN80sim7.dat';
figurename='SG_kh_h5rho_snapshots_old.pdf';

%file= 'basegoryksrhetN80sim2.dat';
%figurename='SG_kc_h5rho_snapshots_old.pdf';



tpoints=[10 60 120 300]/10;




%full figure position
Le=200;Bo=100;Wi=700;He=Wi/4+20;
%[loffset+w*(1+fracsep)*(i-1) b fracw/ncol h] % left, bottom, width, height.
loffset=0.02;
fracsep=0.1;
b=0.05; 
fracw=0.9;
ncol=length(tpoints);
w=fracw/ncol;
%for a squared subplot: h*H=W*w
h=w*Wi/He;

%close all;

tseries1l=load(strcat(folder,file));

N=sqrt(size(tseries1l,2));
nspecies=1;
Nf=size(tseries1l,1)/nspecies;
L=8;
dx=L/N;
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
centroids=[];

h_fig=figure(2);

set(h_fig,'position',[Le,Bo,Wi,He]);
for i = 1:length(tpoints)
    subplot(1,length(tpoints),i)        
    tot42T = squeeze(tseries(1,tpoints(i),:,:));	
    scatter(xst(:,tpoints(i)),yst(:,tpoints(i)),3,'red','filled')
    box on;
    set(gca,'xtick',[],'ytick',[])
    %set(gca,'xticklabel',[])
    %surf(X,Y,-tot42T/dx^2,'FaceColor','interp','EdgeColor', 'none');    
    %colormap gray
    %max(max(tot42T/dx^2))
    xlim([0 L])
    ylim([0 L])
    %zlim([0 zmax]);
    %view(0,90);
    %set(gca,'XColor', 'none','YColor','none')
    %caxis([-2500 0])
    title({[num2str(tpoints(i)*10) ' s']},'FontSize',14)
    %colorbar
    h_subplot(i) = subplot(1,length(tpoints),i);
    %set(h_subplot(i), 'position', [loffset+w*(1+fracsep)*(i-1) b w h] );
    %get(h_subplot(i+(j-1)*length(tpoints)), 'position' )
end
for i = 1:length(tpoints)    
    set(h_subplot(i), 'position', [loffset+w*(1+fracsep)*(i-1) b w h] );
end


saveas(h_fig,figurename)
%saveas(h_fig,'snapshots_old.png')