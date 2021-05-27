
%THIS IS FIGURE 2B
folder="./basegoryrho20nm_gef/"
%Load centered distributions (from center_patch.py)
file= 'output_pos_gef700n1.txt';
%sampling from smoldyn should be every 10s.
tpoints=[10 60 120 300]/10;



Nf=30;
L=8;

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

maxpart = 5000;

timeanalysis=[1:Nf];

xst=-ones(maxpart,length(timeanalysis));
yst=-ones(maxpart,length(timeanalysis));
positions = read_molPos(strcat(folder,file), 30);

line=1;
for i=2:size(positions.Cdc42T,1) %first timpoint (0) is repated, 
    xs=[];
    ys=[];
    
    if size(positions.Cdc42T{i},1) > 0 
        xs =[positions.Cdc42T{i}(:,1) ];
        ys =[positions.Cdc42T{i}(:,2) ];
    end
    if size(positions.BemGEF42{i},1) > 0 
        xs =[xs ; positions.BemGEF42{i}(:,1)];
        ys =[ys ; positions.BemGEF42{i}(:,2)];
    end
    if length(xs) >0
        xst(1:length(xs),line)=xs;
        yst(1:length(ys),line)=ys;
    end 
    line=line+1;

end

xst(xst==-1)=nan;
yst(yst==-1)=nan;	

%TIMESERIES OF PROTEIN DISTRIBUTION
h_fig=figure(1);
%set(gca,'DefaultAxesFontSize',24)

set(h_fig,'position',[Le,Bo,Wi,He]);
for i = 1:length(tpoints)
    subplot(1,length(tpoints),i)        
    scatter(xst(:,tpoints(i)),yst(:,tpoints(i)),3,'red','filled')
    box on;
    set(gca,'xtick',[],'ytick',[])
    xlim([0 L])
    ylim([0 L])
   
    title({[num2str(tpoints(i)*10) ' s']},'FontSize',14)
    %colorbar
    h_subplot(i) = subplot(1,length(tpoints),i);
    %set(h_subplot(i), 'position', [loffset+w*(1+fracsep)*(i-1) b w h] );
    %get(h_subplot(i+(j-1)*length(tpoints)), 'position' )
end
for i = 1:length(tpoints)    
    set(h_subplot(i), 'position', [loffset+w*(1+fracsep)*(i-1) b w h] );
end


saveas(h_fig,'PB_snapshots_old_SI_bifurcation.pdf')
%saveas(h_fig,'snapshots_old.png')
