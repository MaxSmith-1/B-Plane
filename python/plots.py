import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from pathlib import Path
import math


output_directory = Path("output/trials")
names = []

# Sun sphere
r_sun = 695700 
points = np.linspace(0, 2*np.pi, 1000)
x_sun = r_sun*np.cos(points)
y_sun = r_sun*np.sin(points)


plt.style.use('dark_background')


for i in range(1, 4):
    plt.figure(i)
    plt.plot(x_sun, y_sun, color='yellow')


print(output_directory.iterdir())

for i, file in enumerate(output_directory.iterdir()):

    print(type(file))
    
    try:
        names.append(str(file).split("\\")[2].split("_output.csv")[0])
        print(names)
        sim_df = pd.read_csv(file)
    except(PermissionError):
        continue 

    

    print(sim_df)
    


    # ICRF X-Y
    plt.figure(1)
    plt.title("ICRF X vs. Y plane")    
    plt.plot(sim_df.ICRF_X, sim_df.ICRF_Y)
    plt.xlabel("X [km]")
    plt.ylabel("Y [km]")
    plt.axis('square')
    plt.legend("Earth")
    plt.legend(["Sun"] + names)

    # ICRF X-Z
    plt.figure(2)
    plt.title("ICRF X vs. Z plane")
    plt.plot(sim_df.ICRF_X, sim_df.ICRF_Z)
    plt.xlabel("X [km]")
    plt.ylabel("Z [km]")
    plt.axis('square')
    plt.legend(["Sun"] + names)

    # ICRF y-Z
    plt.figure(3)
    plt.title("ICRF Y vs. Z plane")    
    plt.plot(sim_df.ICRF_Y, sim_df.ICRF_Z)
    plt.xlabel("y [km]")
    plt.ylabel("Z [km]")
    plt.axis('square')
    plt.legend(["Sun"] + names)
    

    # Classical Orbital Elements vs. time 
    plt.figure(4)
    plt.title("Semi-major axis vs. time")
    plt.plot(sim_df.time, sim_df.a)
    plt.xlabel("time [s]")
    plt.ylabel("a [km]")
    plt.legend(names)

    plt.figure(5)
    plt.title("Eccentricity vs. time")
    plt.plot(sim_df.time, sim_df.e)
    plt.xlabel("time [s]")
    plt.ylabel("e [-]")
    plt.legend(names)

    plt.figure(6)
    plt.title("Inclination vs. time")
    plt.plot(sim_df.time, sim_df.i)
    plt.xlabel("time [s]")
    plt.ylabel("i [rad]")
    plt.legend(names)

    plt.figure(7)
    plt.title("Longitude of Ascending Node vs. time")
    plt.plot(sim_df.time, sim_df.laan)
    plt.xlabel("time [s]")
    plt.ylabel("Ω [rad]")
    plt.legend(names)

    plt.figure(8)
    plt.title("Argument of Periapsis vs. time")
    plt.plot(sim_df.time, sim_df.gamma)
    plt.xlabel("time [s]")
    plt.ylabel("ω [rad]")
    plt.legend(names)

    plt.figure(9)
    plt.title("True anomaly vs. time")
    plt.plot(sim_df.time, sim_df.f)
    plt.xlabel("time [s]")
    plt.ylabel("f [rad]")
    plt.legend(names)

    plt.figure(10)

    sim_df["v_mag"] = np.sqrt(sim_df["Vx"] ** 2 + sim_df["Vy"] ** 2 + sim_df["Vz"] ** 2)
    plt.title("Velocity magnitude vs. time")
    plt.plot(sim_df.time, sim_df.v_mag)
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [km/s]")
    plt.legend(names)

    plt.figure(11)
    plt.title("Vx vs. time")
    plt.plot(sim_df.time, sim_df.Vx)
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [km/s]")
    plt.legend(names)

    plt.figure(12)
    plt.title("Vy vs. time")
    plt.plot(sim_df.time, sim_df.Vy)
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [km/s]")
    plt.legend(names)

    plt.figure(13)
    plt.title("Vz vs. time")
    plt.plot(sim_df.time, sim_df.Vz)
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [km/s]")
    plt.legend(names)



ref_df = pd.read_csv("python/ref_traj.csv")

target_df = pd.read_csv("output/planets/499.csv")
print(target_df.index)
target_df.columns = ["empty", "date", "X", "Y", "Z", "Vx", "Vy", "Vz", "emprt"]


print(target_df)

plt.figure(1)
plt.plot(ref_df.X, ref_df.Y)
plt.scatter(target_df.X, target_df.Y, c='red')
plt.legend(["Sun"] + names + ["reference trajectory", "Mars"])
ref_df["v_mag"] = np.sqrt(ref_df["VX"] ** 2 + ref_df["VY"] ** 2 + ref_df["VZ"] ** 2)

plt.figure(2)
plt.plot(ref_df.X, ref_df.Z)
plt.scatter(target_df.X, target_df.Z, c='red')
plt.legend(["Sun"] + names + ["reference trajectory", "Mars"])

plt.figure(3)
plt.plot(ref_df.Y, ref_df.Z)
plt.scatter(target_df.Y, target_df.Z, c='red')
plt.legend(["Sun"] + names + ["reference trajectory", "Mars"])




plt.figure(10)
plt.plot(ref_df.time, ref_df.v_mag)
plt.legend(names + ["reference trajectory"])

plt.figure(11)
plt.plot(ref_df.time, ref_df.VX)
plt.legend(names + ["reference trajectory"])

plt.figure(12)
plt.plot(ref_df.time, ref_df.VY)
plt.legend(names + ["reference trajectory"])


plt.figure(13)
plt.plot(ref_df.time, ref_df.VZ)
plt.legend(names + ["reference trajectory"])



plt.show()


