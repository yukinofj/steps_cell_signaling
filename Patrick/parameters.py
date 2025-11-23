import pandas as pd

p = {
    "time step" : 0.1, # 0.25
    "endtime" : 5, # 1
    "DC" : 4e-12,
    "k[0]" : 1e8,
    "k[666]" : 0.3,
    "bc": 4.5,
    # initial conditions large model
    #"big_model_mini_sph_df_path": "/home/yukinofj/code/steps_cell_signaling/Patrick/data_big_model_mini_sph.xls"
    "big_model_mini_sph_df_path": "Patrick/data_big_model_mini_sph.xls"
}

