import pandas as pd
from parameters import p
import numpy as np



val = np.linspace(0, 2, 10)

def format_param(val, decimals=4):
    return float(f"{val:.{decimals}g}")

for idx, value in enumerate(val):
    print(idx, format_param(value))

#print(list(format_param(val)))