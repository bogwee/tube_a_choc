import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("output.csv")
x = df['x'].to_numpy()
rho = df['rho'].to_numpy()
U = df['U'].to_numpy()
p = df['p'].to_numpy()

plt.figure(figsize=(15, 4))

plt.subplot(1, 3, 1)
plt.plot(x, rho, label='Density (rho)', color='blue')
plt.title('Masse Volumique (t=0.2s)')
plt.xlabel('Position (x)')
plt.ylabel('Masse Volumique (rho)')

plt.subplot(1, 3, 2)
plt.plot(x, U, label='Velocity (U)', color='orange')
plt.title('Vitesse (t=0.2s)')
plt.xlabel('Position (x)')
plt.ylabel('Vitesse (U)')

plt.subplot(1, 3, 3)
plt.plot(x, p, label='Pressure (P)', color='green')
plt.title('Pression (t=0.2s)')
plt.xlabel('Position (x)')
plt.ylabel('Pression (P)')

plt.show()

