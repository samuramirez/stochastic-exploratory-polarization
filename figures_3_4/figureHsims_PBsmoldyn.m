function figureHsims_PBsmoldyn(gefconc) 
%6 concentrations, 3 methods, 3 grid sizes, 30sims, 300s-10sinterval 30 timepoints, 50 Hvals
%~22MB

%Check in function computing H
%max_r  = 5;
%dr=0.1;
Nsims=30;
%GEFconc=["200","300","400","500","600","700"];
%GEFconc=["500","600","700"];
GEFconc=num2str(gefconc)
%timeanalysis=[1:Nf];%verify this with H function file
    Hallsims=[];
    GEFconc
    for l = 1: Nsims
    %for l = 3    
        l
        filename=strcat('./basegoryrho20nm_gef','/output_pos_gef',GEFconc,'n',num2str(l),'.txt');
        [H_vals,r_vals] = Hfn_smoldyn_clustQuant(filename); %can we compute this just for r = 0.5 ? 
        %Hsr0p5 = [Hsr0p5 ; H_vals(find(r_vals==0.5),:)];
        Hallsims(l,:,:)=H_vals; %nsims x 50Hvalues * 30timepoints 
    end
    save(strcat('./HvaluesPB/Hvalues','gef',GEFconc, 'PB.mat'), 'Hallsims');
 

%figure(1)
%hold on;
%for j = 1: Nsims
%	plot(Hsr0p5(j,:))
%end
%hold off;
end
