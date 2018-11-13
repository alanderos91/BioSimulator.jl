import stochpy
import statistics as stats
from statsmodels import robust
import argparse
import random
import sys

def benchmark(fdir, method, tfinal, nsaves, nreal, seed, n_sample):
  # initialize
  smod = stochpy.SSA()
  times = []

  smod.Model(model_file='stochpy.psc', dir=fdir)

  random.seed(seed)

  for i in range(n_sample):
    print(f"sample {i+1} / {n_sample}")
    smod.DoStochSim(trajectories=nreal, end=tfinal, method=method, mode="time")
    times.append(smod.simulation_time)

  return times

def save_results(times, model_name):
  fout = open(f"./stochpy/{model_name}.txt", 'a')
  
  for t in times:
    fout.write("%s\n" % t)

  fout.close()

  return

parser = argparse.ArgumentParser()

parser.add_argument(dest='model_name')
parser.add_argument(dest='tfinal', type=float)
parser.add_argument(dest='nsaves', type=int)
parser.add_argument(dest='nreal', type=int)
parser.add_argument(dest='seed', type=int)
parser.add_argument(dest='n_sample', type=int)

args = parser.parse_args()

# METHODS = ['Direct', 'FRM', 'NRM', 'TauLeap']
METHODS = ['Direct']

fdir     = f"./{args.model_name}/"
tfinal   = args.tfinal
nsaves   = args.nsaves
nreal    = args.nreal
seed     = args.seed
n_sample = args.n_sample

for method in METHODS:
  times = benchmark(fdir, method, tfinal, nsaves, nreal, seed, n_sample)

  save_results(times, args.model_name)
