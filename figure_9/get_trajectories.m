%Get the tracjectories of the center of mass of active Cdc42 projected on the sphere

folder='D:/dynamic_polarity_data/fig_3D_transition_quant/config_files_quant/';
nameroot = '3d_newreac_gef100_k4a_';
params=["0p001","0p002","0p005","0p01","0p02","0p05","0p1","1"];
% name example: 3d_newreac_gef100_k4a_0p001_01.xyz
extension_in='.xyz';
extension_out='.txt';
Nsims=30;
maxframes=35; %every frame = 1 min, most simulations didn't get to 20 min
R=4.5135/2;

for param = params
    for i = 1:Nsims
        realization=sprintf('%02d',i);
        filenamebase=[folder nameroot char(param) '_' realization];
        [t,trajectory]=calculate_activeCOM_trajectory([filenamebase extension_in ],maxframes,R);
        filenameout =[folder 'traj_' nameroot char(param) '_' realization extension_out];
        save(filenameout,'trajectory', '-ascii');
    end
end

