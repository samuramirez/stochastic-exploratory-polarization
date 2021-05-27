% 2020 May 12 / Mike Pablo
% Generates a smoldyn configuration file for 3D particle-based simulations of
% yeast polarity establishment, using Sam Ramirez's biochemical network.
% ---------------------------------------------------------------------------
% INPUTS
% fileprefix    - String. Used to name the configuration file.
% tstop         - Double. Desired final simulation time (seconds).
% dtime        - Double. Time step (seconds).
% samplingrate  - Double. How often, in numbers of timesteps, to record data
% random_seed   - Sets random seed in smoldyn, enables reproducibility.
% k1a*          - First order rate constant (s-1): BemGEF(c) -> BemGEF(m)
% k1b           - First order rate constant (s-1): BemGEF(m) -> BemGEF(c)
% P2a           - Second order probability rate (prob/timestep): Cdc42D(m) + BemGEF(m) -> Cdc42T(m) + BemGEF(m)
% k2b           - First order rate constant (s-1): Cdc42T(m) -> Cdc42D(m)
% P3            - Second order probability rate (prob/timestep): Cdc42D(m) + BemGEF42(m) -> Cdc42T(m) + BemGEF42(m)
% P4a           - Second order probability rate (prob/timestep): Cdc42T(m) + BemGEF(m)   -> BemGEF42(m)
% k4b           - First order rate constant (s-1): BemGEF42(m) -> BemGEF(m) + Cdc42T(m)
% k5a*          - First order rate constant (s-1): Cdc42D(c) -> Cdc42D(m) [PM association rate, see note]
% k5b           - First order rate constant (s-1): Cdc42D(m) -> Cdc42D(c)
% P7            - Second order probability rate (prob/timestep): Cdc42T(m) + BemGEF(c)   -> BemGEF42(m)
% P8            - Second order probability rate (prob/timestep): Cdc42D(c) + BemGEF42(m) -> Cdc42T(m) + BemGEF42(m)
% Dm, Dc        - Diffusion coefficients of membrane and cytosol, respectively (um2/s)
% rho           - Reaction radius for particles (um)
% d_sphere      - Diameter of the simulation volume (um)
% n_BemGEF      - Number of BemGEF molecules
% n_Cdc42       - Number of Cdc42 molecules
%
% * [PM association rate, see a note below in the code]
% 
%
% OUTPUT
% A .cfg file is written to the current folder with name fileprefix.cfg
% Running this .cfg file with smoldyn will generate data with names
% fileprefix.xyz, containing molecular coordinates over time, and
% fileprefix.num. See analysis code for more details.
% ------------------------------------------------------------------

function make_smoldyn_cfg(fileprefix,tstop,dtime,samplingrate,random_seed,...
                          k1a,k1b,P2a,k2b,P3,P4a,k4b,k5a,k5b,P7,P8,...
                          Dm,Dc,rho,d_sphere,n_BemGEF,n_Cdc42)

cfg_name = [fileprefix '.cfg']; % configuration file for smoldyn
xyz_name = [fileprefix '.xyz']; % xyz coords of species
num_name = [fileprefix '.num']; % number of species

fid=fopen(cfg_name,'w');

% ---------------------------------------------------------------------------
% PM association rates, in 3D smoldyn simulations, must be specified as
% "adsorption coefficients" with units of length / time. These coefficients
% can be computed from a PM association rate by:
%     adsorption coeff = rate * z,
% where z = V_simulation_domain / Area_simulation_domain
%       --> if d_sphere = 4.5135 um, z = 0.7523 um.
% The connection to the first-order rate in 2D simulations breaks down in
% the diffusion-limited context (i.e. a poorly mixed cytoplasm).
% ----------------------------------------------------------------------------
z = (d_sphere/2)/3; % Ratio between sphere volume and sphere surface area.

fprintf(fid,'random_seed %i\n',random_seed);

fprintf(fid,'variable d_sphere = %g\n',d_sphere);
fprintf(fid,'variable r_sphere = d_sphere/2\n');
fprintf(fid,'variable d_sphereplus_1Leps = d_sphere+0.1\n');
fprintf(fid,'variable d_sphereplus_2Leps = d_sphere+2*0.1\n');
fprintf(fid,'variable rho = %g\n',rho);
fprintf(fid,'variable rho_eps = 0.00001\n');

fprintf(fid,'variable k1a = %g\n',k1a);
fprintf(fid,'variable k1b = %g\n',k1b);
fprintf(fid,'variable P2a = %g\n',P2a);
fprintf(fid,'variable k2b = %g\n',k2b);
fprintf(fid,'variable P3 = %g\n',P3);
fprintf(fid,'variable P4a = %g\n',P4a);
fprintf(fid,'variable k4b = %g\n',k4b);
fprintf(fid,'variable k5a = %g\n',k5a);
fprintf(fid,'variable k5b = %g\n',k5b);
fprintf(fid,'variable P7 = %g\n',P7);
fprintf(fid,'variable P8 = %g\n',P8);

fprintf(fid,'dim 3\n');
fprintf(fid,'boundaries x -0.1 d_sphereplus_1Leps\n');
fprintf(fid,'boundaries y -0.1 d_sphereplus_1Leps\n');
fprintf(fid,'boundaries z -0.1 d_sphereplus_1Leps\n');

fprintf(fid,'species Cdc42T Cdc42D BemGEF BemGEF42 complex_Cdc42Dm_BemGEF42 complex_Cdc42Dm_BemGEFm\n');
fprintf(fid,'difc Cdc42T(all) %g\n',Dm);
fprintf(fid,'difc Cdc42D(soln) %g\n',Dc);
fprintf(fid,'difc Cdc42D(down) %g\n',Dm);
fprintf(fid,'difc BemGEF42(all) %g\n',Dm);
fprintf(fid,'difc BemGEF(soln) %g\n',Dc);
fprintf(fid,'difc BemGEF(down) %g\n',Dm);
fprintf(fid,'difc complex_Cdc42Dm_BemGEFm(all) %g\n',Dm);
fprintf(fid,'difc complex_Cdc42Dm_BemGEF42(all) %g\n',Dm);

fprintf(fid,'molecule_lists list1 list2 list3 list4\n');
fprintf(fid,'mol_list Cdc42T(all) list1\n');
fprintf(fid,'mol_list Cdc42D(all) list2\n');
fprintf(fid,'mol_list BemGEF42(all) list3\n');
fprintf(fid,'mol_list BemGEF(all) list4\n');
fprintf(fid,'variable z = %g\n',z);

fprintf(fid,'start_surface inner_walls\n');
fprintf(fid,'action both all reflect\n');
fprintf(fid,'polygon both edge\n');
fprintf(fid,'rate BemGEF soln down k1a*z\n');
fprintf(fid,'rate BemGEF down soln k1b\n');
fprintf(fid,'rate Cdc42D soln down k5a*z\n');
fprintf(fid,'rate Cdc42D down soln k5b\n');
fprintf(fid,'panel sphere r_sphere r_sphere r_sphere -r_sphere 50 50\n');
fprintf(fid,'end_surface\n');

fprintf(fid,'start_surface outer_walls\n');
fprintf(fid,'action both all reflect\n');
fprintf(fid,'polygon both none\n');
fprintf(fid,'panel rect +x -0.1 -0.1 -0.1 d_sphereplus_2Leps d_sphereplus_2Leps\n');
fprintf(fid,'panel rect -x d_sphereplus_2Leps -0.1 -0.1 d_sphereplus_2Leps d_sphereplus_2Leps\n');
fprintf(fid,'panel rect +y -0.1 -0.1 -0.1 d_sphereplus_2Leps d_sphereplus_2Leps\n');
fprintf(fid,'panel rect -y -0.1 d_sphereplus_2Leps -0.1 d_sphereplus_2Leps d_sphereplus_2Leps\n');
fprintf(fid,'panel rect +z -0.1 -0.1 -0.1 d_sphereplus_2Leps d_sphereplus_2Leps\n');
fprintf(fid,'panel rect -z -0.1 -0.1 d_sphereplus_2Leps d_sphereplus_2Leps d_sphereplus_2Leps\n');
fprintf(fid,'end_surface\n');

fprintf(fid,'start_compartment full_domain\n');
fprintf(fid,'surface inner_walls\n');
fprintf(fid,'point r_sphere r_sphere r_sphere\n');
fprintf(fid,'end_compartment\n');

fprintf(fid,'reaction Cdc42Dm_2_T_bindTo_BemGEFm  Cdc42D(down) + BemGEF(down) -> complex_Cdc42Dm_BemGEFm(down)\n');
fprintf(fid,'reaction Cdc42Dm_2_T_catBy_BemGEFm complex_Cdc42Dm_BemGEFm(down) -> Cdc42T(down) + BemGEF(down)\n');
fprintf(fid,'reaction_probability Cdc42Dm_2_T_bindTo_BemGEFm P2a\n');
fprintf(fid,'binding_radius Cdc42Dm_2_T_bindTo_BemGEFm rho\n');
fprintf(fid,'reaction_probability Cdc42Dm_2_T_catBy_BemGEFm 1\n');
fprintf(fid,'product_placement Cdc42Dm_2_T_catBy_BemGEFm unbindrad rho+rho_eps # place just beyond region\n');

fprintf(fid,'reaction Cdc42T_2_Cdc42Dm Cdc42T(down) -> Cdc42D(down) k2b\n');

fprintf(fid,'reaction Cdc42Dm_2_T_bindTo_BemGEF42  Cdc42D(down) + BemGEF42(down) -> complex_Cdc42Dm_BemGEF42(down)\n');
fprintf(fid,'reaction Cdc42Dm_2_T_catBy_BemGEF42 complex_Cdc42Dm_BemGEF42(down) -> Cdc42T(down) + BemGEF42(down)\n');
fprintf(fid,'reaction_probability Cdc42Dm_2_T_bindTo_BemGEF42 P3\n');
fprintf(fid,'binding_radius Cdc42Dm_2_T_bindTo_BemGEF42 rho\n');
fprintf(fid,'reaction_probability Cdc42Dm_2_T_catBy_BemGEF42 1\n');
fprintf(fid,'product_placement Cdc42Dm_2_T_catBy_BemGEF42 unbindrad rho+rho_eps # place just beyond region\n');

fprintf(fid,'reaction Cdc42Dc_2_T_bindTo_BemGEF42  Cdc42D(soln) + BemGEF42(down) -> complex_Cdc42Dm_BemGEF42(down)\n');
fprintf(fid,'reaction_probability Cdc42Dc_2_T_bindTo_BemGEF42 P8\n');
fprintf(fid,'binding_radius Cdc42Dc_2_T_bindTo_BemGEF42 rho\n');

fprintf(fid,'reaction make_BemGEF42_fromm BemGEF(down) + Cdc42T(down) <-> BemGEF42(down)\n');
fprintf(fid,'reaction_probability make_BemGEF42_frommfwd P4a\n');
fprintf(fid,'binding_radius make_BemGEF42_frommfwd rho\n');

fprintf(fid,'reaction_rate make_BemGEF42_frommrev k4b\n');
fprintf(fid,'product_placement make_BemGEF42_frommrev unbindrad rho+rho_eps\n');

fprintf(fid,'reaction make_BemGEF42_fromc BemGEF(soln) + Cdc42T(down) -> BemGEF42(down)\n');
fprintf(fid,'reaction_probability make_BemGEF42_fromc P7\n');
fprintf(fid,'binding_radius make_BemGEF42_fromc rho\n');

fprintf(fid,'time_start 0\n');
fprintf(fid,'time_stop %i\n',tstop);
fprintf(fid,'time_step %g\n',dtime);

fprintf(fid,'compartment_mol %i Cdc42D(soln) full_domain\n',n_Cdc42);
fprintf(fid,'compartment_mol %i BemGEF(soln) full_domain\n',n_BemGEF);

fprintf(fid,'output_files %s\n',num_name);
fprintf(fid,'cmd B molcountheader %s\n',num_name);
fprintf(fid,'cmd N %i molcount %s\n',samplingrate,num_name);

fprintf(fid,'output_files %s\n',xyz_name);
fprintf(fid,'cmd B molpos Cdc42T(all) %s\n',xyz_name);
fprintf(fid,'cmd B molpos BemGEF42(all) %s\n',xyz_name);
fprintf(fid,'cmd B molpos Cdc42D(down) %s\n',xyz_name);
fprintf(fid,'cmd B molpos Cdc42D(soln) %s\n',xyz_name);
fprintf(fid,'cmd B molpos BemGEF(down) %s\n',xyz_name);
fprintf(fid,'cmd B molpos BemGEF(soln) %s\n',xyz_name);
fprintf(fid,'cmd B molpos complex_Cdc42Dm_BemGEFm(all) %s\n',xyz_name);
fprintf(fid,'cmd B molpos complex_Cdc42Dm_BemGEF42(all) %s\n',xyz_name);
fprintf(fid,'cmd N %i molpos Cdc42T(all) %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos BemGEF42(all) %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos Cdc42D(down) %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos Cdc42D(soln) %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos BemGEF(down) %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos BemGEF(soln) %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos complex_Cdc42Dm_BemGEFm(all) %s\n',samplingrate,xyz_name);
fprintf(fid,'cmd N %i molpos complex_Cdc42Dm_BemGEF42(all) %s\n',samplingrate,xyz_name);
fclose(fid);
end
