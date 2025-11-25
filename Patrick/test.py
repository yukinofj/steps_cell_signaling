import pandas as pd
from parameters import p
import numpy as np
import argparse

pars = argparse.ArgumentParser()
pars.add_argument("runnr", type=int, help="Run number")
pars.add_argument("n_reps", type=int, help="Number of repetitions")
pars.add_argument("ellipsoidity", type=float, help="Ellipsoidity")

args = pars.parse_args()

print("runnr:", args.runnr)
print("n_reps:", args.n_reps)
print("ellipsoidity:", args.ellipsoidity)

runnr, n_reps, ellipsoidity = (args.runnr, args.n_reps, args.ellipsoidity)

print(runnr, n_reps, ellipsoidity)