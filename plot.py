import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

output = pd.read_csv('output.csv')
rho = output['rho'].values.to_numpy()
U = output['U'].values.to_numpy()
p = output['P'].values.to_numpy()
np.array(output)

