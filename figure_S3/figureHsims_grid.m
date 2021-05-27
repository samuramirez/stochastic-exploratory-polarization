function  figureHsims_grid(folder,variable,param,Nsims)

%6 concentrations, 3 methods, 3 grid sizes, 30sims, 300s-10sinterval 30 timepoints, 50 Hvals
%~22MB

%Check in function computing H
%max_r  = 5;
%dr=0.1;


%display(meth)
%timeanalysis=[1:Nf];%verify this with H function file
Hallsims=[];
for l = 1: Nsims
%for l = 1    
    l
    filename=strcat(folder,'/dyn','N80',variable,param,'sim',num2str(l),'.dat');
    [H_vals,r_vals] = Hfn_grid_clustQuant(filename); %can we compute this just for r = 0.5 ? 
    %Hsr0p5 = [Hsr0p5 ; H_vals(find(r_vals==0.5),:)];
    Hallsims(l,:,:)=H_vals; %nsims x 50Hvalues * 30timepoints 
end
save(strcat(folder,'/Hvalues',variable,param,'N80','.mat'), 'Hallsims');

%figure(1)
%hold on;
%for j = 1: Nsims
%	plot(Hsr0p5(j,:))
%end
%hold off;

end
