import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import math

# Plot b-plane here
output_directory = Path("output")
names = []

mars_r = 514750 # km

for i, file in enumerate(output_directory.iterdir()):

    print(type(file))
    
    try:
        names.append(str(file).split("\\")[1].split("_output.csv")[0])
        sim_df = pd.read_csv(file)
    except(PermissionError):
        continue 

    

    print(sim_df)

    plt.figure(1)
    colors = np.arange(len(sim_df["b_impact_parameter_x"]))

    plt.scatter(sim_df["b_impact_parameter_x"], sim_df["b_impact_parameter_x"], c=colors, cmap="rainbow")

    plt.figure(2)
    plt.scatter(sim_df["time"], sim_df["V_infinity"])

    plt.figure(3)
    plt.scatter(sim_df["time"], sim_df["r_soi"])

    plt.figure(4)
    plt.scatter(sim_df["time"], sim_df["b_impact_parameter"])

    plt.scatter(sim_df["time"], sim_df["r_soi"])




mars_x = mars_r * np.cos(np.linspace(0, 2*np.pi, 100))
mars_y = mars_r * np.sin(np.linspace(0, 2*np.pi, 100))

plt.figure(1)
plt.plot(mars_x, mars_y, 'red')
plt.title("B-Plane")
plt.xlabel("B_X [km]")
plt.ylabel("B_Y [km]")

plt.figure(2)
plt.title("V infinity")
plt.ylabel("V_infinity [km / s]")
plt.xlabel("time [s]")


plt.show()
