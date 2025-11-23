import pandas as pd
from parameters import p


fraction = [0.0, 0.5, 1.0]
base_count = 1000
for c_idx, fraction in enumerate(fraction):
    c_scaled = fraction * base_count
    print(c_scaled)