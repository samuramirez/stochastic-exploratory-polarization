function  figureHsims_grid(gefconc,nmeth)

%6 concentrations, 3 methods, 3 grid sizes, 30sims, 300s-10sinterval 30 timepoints, 50 Hvals
%~22MB

%Check in function computing H
%max_r  = 5;
%dr=0.1;
Nsims=30;
%meth=['ka','ksrhet','kh'];
%{
if ktype==4
	meth='ksrhet';
	display(meth)
	display(gefconc)
	display(nmeth)
else if ktype==1
	meth='kh';
	display(meth)
else
	display("wrong ktype")
end
%}

meth='ksrhet';


display(meth)
Nmeth=num2str(nmeth)
GEFconc=num2str(gefconc);
%Nmeth=["40","80","120","160"];
%Nmeth=["80"];
%GEFconc=["200","300","400","500","600","700"];
%GEFconc=["700"];
%timeanalysis=[1:Nf];%verify this with H function file
            Hallsims=[];
            display(GEFconc)
            display(Nmeth)
            for l = 1: Nsims
            %for l = 1    
                l
                filename=strcat('./basegorymeth_rho20nmgef',GEFconc,'kdmeso/basegory',meth,'N',Nmeth,'sim',num2str(l),'.dat');
                [H_vals,r_vals] = Hfn_grid_clustQuant(filename); %can we compute this just for r = 0.5 ? 
                %Hsr0p5 = [Hsr0p5 ; H_vals(find(r_vals==0.5),:)];
                Hallsims(l,:,:)=H_vals; %nsims x 50Hvalues * 30timepoints 
            end
            save(strcat('./HvaluesSG/Hvalues','gef',GEFconc,meth,'N',Nmeth,'.mat'), 'Hallsims');

%figure(1)
%hold on;
%for j = 1: Nsims
%	plot(Hsr0p5(j,:))
%end
%hold off;

end
