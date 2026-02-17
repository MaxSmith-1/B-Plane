import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import math

# Plot b-plane here
output_directory = Path("output/trials")
names = []

mars_r = 3396  # km
plt.style.use('dark_background')
bs_list = []
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
    plt.scatter(sim_df["b_impact_parameter_x"], sim_df["b_impact_parameter_y"], c=colors, cmap="rainbow")

    plt.figure(2)
    plt.scatter(sim_df["time"], sim_df["V_infinity"])

    plt.figure(3)
    plt.scatter(sim_df["time"], sim_df["r_soi"])

    plt.figure(4)
    plt.scatter(sim_df["time"], sim_df["b_impact_parameter"])
    plt.title("Impact parameter magnitude")
    plt.xlabel("time [s]")
    plt.ylabel("Impact Parameter [km]")


    
    b_x_y = list(sim_df[sim_df["passed_b_plane"] == True][["b_impact_parameter_x", "b_impact_parameter_y"]].iloc[0])

    bs_list.append(b_x_y)

    



bs = pd.DataFrame(bs_list, columns=["b_impact_parameter_x", "b_impact_parameter_y"])

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

plt.figure(5)

print(bs)
plt.scatter(bs["b_impact_parameter_x"], bs["b_impact_parameter_y"])
plt.title("Impact Parameter")
plt.xlabel("B_X [km]")
plt.ylabel("B_Y [km]")
plt.plot(mars_x, mars_y, 'red')

# Draw covariance sphere




import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path

# ... [Your existing setup code] ...

# Assuming we are animating the last file read in your loop
# Adjust indexing if you want to animate a specific trial
df = sim_df 

fig, ax = plt.subplots(figsize=(8, 6))

# Setup the B-Plane background (Mars/Target)
mars_r = 514750 
theta = np.linspace(0, 2*np.pi, 100)
ax.plot(mars_r * np.cos(theta), mars_r * np.sin(theta), 'red', label='Target SOI')

# Initialize the line and point
line, = ax.plot([], [], 'b-', alpha=0.5)
scat = ax.scatter([], [], c=[], cmap="rainbow", vmin=0, vmax=len(df))

# Set plot limits based on data
ax.set_xlim(df["b_impact_parameter_x"].min() * 1.1, df["b_impact_parameter_x"].max() * 1.1)
ax.set_ylim(df["b_impact_parameter_y"].min() * 1.1, df["b_impact_parameter_y"].max() * 1.1)
ax.set_title("B-Plane Trajectory Evolution")
ax.set_xlabel("B_X [km]")
ax.set_ylabel("B_Y [km]")
ax.grid(True, linestyle='--', alpha=0.6)

def init():
    line.set_data([], [])
    return line,

def update(frame):
    # Get data up to current frame
    x = df["b_impact_parameter_x"].iloc[:frame]
    y = df["b_impact_parameter_y"].iloc[:frame]
    
    line.set_data(x, y)
    
    # Update scatter to show the "head" of the trajectory
    # Using 'c' to maintain your rainbow color evolution
    scat.set_offsets(np.column_stack([x, y]))
    scat.set_array(np.arange(frame))
    
    return line, scat

# Create animation (frames=number of rows, interval=ms between frames)
ani = FuncAnimation(fig, update, frames=len(df), init_func=init, blit=True, interval=20)

plt.show()
