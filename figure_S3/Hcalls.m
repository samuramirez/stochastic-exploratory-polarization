
folder = "../basegory/basegory_k1a_100x_gef";
folder = "../basegory/basegory_k2a_10_gef";
folder = "../basegory/basegory_k2b_0p0625x_gef";
folder = "../basegory/basegory_k3_10_gef";
folder = "../basegory/basegory_k5a_100x_gef";

folder = "D:/dynamic_polarity_data/basegory/basegory_k2b_0p0625x_k5a_10x_gef";
folder = "D:/dynamic_polarity_data/basegory/basegory_k7_100_gef";
folder = "D:/dynamic_polarity_data/basegory_newreac_gef";
folder = "D:/dynamic_polarity_data/basegory/basegory_k2b_0p0625x_k5a_10x_gef_long";



variable = "gef";
Nsims=5;

for param = ["15","25","50"]
    figureHsims_grid(folder,variable,param,Nsims)
    
end