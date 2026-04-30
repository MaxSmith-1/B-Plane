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
    # plt.legend(["Sun"])


target_df = pd.read_csv("output/planets/499.csv")

target_df.columns = ["empty", "date", "X", "Y", "Z", "Vx", "Vy", "Vz", "emprt"]

ref_df = pd.read_csv("python/veo_ref_traj.csv")

for i, file in enumerate(output_directory.iterdir()):


    
    try:
        names.append(str(file).split("\\")[2].split("_output.csv")[0])

        sim_df = pd.read_csv(file)
    except(PermissionError):
        continue 



    # J2000 X-Y
    plt.figure(1)
    plt.title("J2000 X vs. Y plane")    
    plt.scatter(target_df.X, target_df.Y, c='red')
    plt.plot(sim_df.ICRF_X, sim_df.ICRF_Y)
    plt.xlabel("X [km]")
    plt.ylabel("Y [km]")
    plt.axis('square')
    # plt.legend("Earth")
    # plt.legend(["Sun"] + names)
    # plt.legend("Earth")
    # plt.legend(["Sun"] + names)

    # J2000 X-Z
    plt.figure(2)
    plt.title("J2000 X vs. Z plane")
    plt.scatter(target_df.X, target_df.Z, c='red')
    plt.plot(sim_df.ICRF_X, sim_df.ICRF_Z)
    plt.xlabel("X [km]")
    plt.ylabel("Z [km]")
    plt.axis('square')
    # plt.legend(["Sun"] + names)
    # plt.legend(["Sun"] + names)

    # J2000 y-Z
    plt.figure(3)
    plt.title("J2000 Y vs. Z plane")  
    plt.scatter(target_df.Y, target_df.Z, c='red')  
    plt.plot(sim_df.ICRF_Y, sim_df.ICRF_Z, color="white")
    plt.xlabel("Y [km]")
    plt.ylabel("Z [km]")
    plt.axis('square')
    # plt.legend(["Sun"] + names)
    # plt.legend(["Sun"] + names)
    

    # Classical Orbital Elements vs. time 
    plt.figure(4)
    plt.title("Semi-major axis vs. time")
    plt.plot(sim_df.time, sim_df.a)
    plt.xlabel("time [s]")
    plt.ylabel("a [km]")
    # plt.legend(names)
    # plt.legend(names)

    plt.figure(5)
    plt.title("Eccentricity vs. time")
    plt.plot(sim_df.time, sim_df.e)
    plt.xlabel("time [s]")
    plt.ylabel("e [-]")
    # plt.legend(names)
    # plt.legend(names)

    plt.figure(6)
    plt.title("Inclination vs. time")
    plt.plot(sim_df.time, sim_df.i)
    plt.xlabel("time [s]")
    plt.ylabel("i [rad]")
    # plt.legend(names)
    # plt.legend(names)

    plt.figure(7)
    plt.title("Longitude of Ascending Node vs. time")
    plt.plot(sim_df.time, sim_df.laan)
    plt.xlabel("time [s]")
    plt.ylabel("Ω [rad]")
    # plt.legend(names)
    # plt.legend(names)

    plt.figure(8)
    plt.title("Argument of Periapsis vs. time")
    plt.plot(sim_df.time, sim_df.gamma)
    plt.xlabel("time [s]")
    plt.ylabel("ω [rad]")
    # plt.legend(names)
    # plt.legend(names)

    plt.figure(9)
    plt.title("True anomaly vs. time")
    plt.plot(sim_df.time, sim_df.f)
    plt.xlabel("time [s]")
    plt.ylabel("f [rad]")
    # plt.legend(names)
    # plt.legend(names)

    plt.figure(10)

    sim_df["v_mag"] = np.sqrt(sim_df["Vx"] ** 2 + sim_df["Vy"] ** 2 + sim_df["Vz"] ** 2)
    plt.title("Velocity magnitude vs. time")
    plt.plot(sim_df.time, sim_df.v_mag)
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [km/s]")
    # plt.legend(names)
    # plt.legend(names)

    plt.figure(11)
    plt.title("Vx vs. time")
    plt.plot(sim_df.time, sim_df.Vx)
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [km/s]")
    # plt.legend(names)
    # plt.legend(names)

    plt.figure(12)
    plt.title("Vy vs. time")
    plt.plot(sim_df.time, sim_df.Vy)
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [km/s]")
    # plt.legend(names)
    # plt.legend(names)

    plt.figure(13)
    plt.title("Vz vs. time")
    plt.plot(sim_df.time, sim_df.Vz)
    plt.xlabel("Time [s]")
    plt.ylabel("Velocity [km/s]")
    # plt.legend(names)
    # plt.legend(names)






plt.figure(1)

plt.plot(ref_df.X, ref_df.Y, color="#26DFF0")

# plt.legend(["Sun"] + names + ["Mars"])
ref_df["v_mag"] = np.sqrt(ref_df["VX"] ** 2 + ref_df["VY"] ** 2 + ref_df["VZ"] ** 2)

#plt.legend(["Sun", "Mars", "Integrated Trajectory", "Tabulated Trajectory"])

plt.figure(2)

plt.plot(ref_df.X, ref_df.Z,  color="#26DFF0")

# plt.legend(["Sun"] + names + ["Mars"])
#plt.legend(["Sun", "Mars", "Integrated Trajectory", "Tabulated Trajectory"])

plt.figure(3)

plt.plot(ref_df.Y, ref_df.Z,  color="#26DFF0")

# plt.legend(["Sun"] + names + ["Mars"])
#plt.legend(["Sun", "Mars", "Integrated Trajectory", "Tabulated Trajectory"])




plt.figure(10)
plt.plot(ref_df.time, ref_df.v_mag)
#plt.legend(names + ["M.R.O, NASA Horizons Trajectory"])

# plt.figure(11)
# plt.plot(ref_df.time, ref_df.VX)
# plt.legend(names + ["M.R.O, NASA Horizons Trajectory"])

# plt.figure(12)
# plt.plot(ref_df.time, ref_df.VY)
# plt.legend(names + ["M.R.O, NASA Horizons Trajectory"])


# plt.figure(13)
# plt.plot(ref_df.time, ref_df.VZ)
# plt.legend(names + ["M.R.O, NASA Horizons Trajectory"])

x_diff = abs(sim_df["ICRF_X"].iloc[-1] - ref_df["X"].iloc[-1])
print(x_diff)
y_diff = abs(sim_df["ICRF_Y"].iloc[-1] - ref_df["Y"].iloc[-1])
print(y_diff)
z_diff = abs(sim_df["ICRF_Z"].iloc[-1] - ref_df["Z"].iloc[-1])
print(z_diff)

print("Simmed r")
print(np.linalg.norm(np.array([sim_df["ICRF_X"].iloc[-1], sim_df["ICRF_Y"].iloc[-1], sim_df["ICRF_Z"].iloc[-1]])))

print("tabulated r")
print(np.linalg.norm(np.array([ref_df["X"].iloc[-1], ref_df["Y"].iloc[-1], ref_df["Z"].iloc[-1]])))

print("r difference")

print(np.linalg.norm(np.array([x_diff, y_diff, z_diff])))

print("r percent difference")
print(np.linalg.norm(np.array([x_diff, y_diff, z_diff])) / np.linalg.norm(np.array([ref_df["X"].iloc[-1], ref_df["Y"].iloc[-1], ref_df["Z"].iloc[-1]])))


plt.show()

