import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation
from pathlib import Path

# ── Constants ─────────────────────────────────────────────
MARS_R = 3396  # km
SOI_R  = 514_750  # km
COLS   = [          # only load what you actually plot
    "b_impact_parameter_x",
    "b_impact_parameter_y",
    "b_impact_parameter",
    "time",
    "V_infinity",
    "r_soi",
    "passed_b_plane",
]

output_directory = Path("output/trials")
plt.style.use("dark_background")

# ── Pre-allocate accumulators (avoid repeated DataFrame concat) ──
bs_list   = []          # final b-plane hits
all_bx    = []          # figure-1 / figure-4 x data
all_by    = []
all_time  = []
all_vinf  = []
all_rsoi  = []
all_bmag  = []
all_color_offsets = []  # cumulative index for rainbow colour
cumulative = 0

last_df = None           # keep only the last df for animation

mars_circle = np.exp(1j * np.linspace(0, 2 * np.pi, 100)) * MARS_R  # complex trick → 1 alloc

# ── Main loop ─────────────────────────────────────────────
for file in output_directory.iterdir():
    try:
        sim_df = pd.read_csv(file, usecols=COLS)   # ← only needed columns
    except (PermissionError, ValueError):
        continue

    n = len(sim_df)

    # Accumulate numpy arrays instead of re-drawing every iteration
    all_bx.append(sim_df["b_impact_parameter_x"].to_numpy())
    all_by.append(sim_df["b_impact_parameter_y"].to_numpy())
    all_time.append(sim_df["time"].to_numpy())
    all_vinf.append(sim_df["V_infinity"].to_numpy())
    all_rsoi.append(sim_df["r_soi"].to_numpy())
    all_bmag.append(sim_df["b_impact_parameter"].to_numpy())
    all_color_offsets.append(np.arange(cumulative, cumulative + n))
    cumulative += n

    hit = sim_df[sim_df["passed_b_plane"].astype(bool)]
    if not hit.empty:
        bs_list.append(hit[["b_impact_parameter_x", "b_impact_parameter_y"]].iloc[0].tolist())

    last_df = sim_df   # only reference, old df gets GC'd next iteration

# ── Concatenate once (fast, single allocation) ────────────
bx      = np.concatenate(all_bx)
by      = np.concatenate(all_by)
time_   = np.concatenate(all_time)
vinf    = np.concatenate(all_vinf)
rsoi    = np.concatenate(all_rsoi)
colors  = np.concatenate(all_color_offsets)
bmag = np.concatenate(all_bmag)


print(len(all_bx))

# Free the per-file lists immediately
del all_bx, all_by, all_time, all_vinf, all_rsoi, all_color_offsets, all_bmag

# ── Static figures (draw once, not N times) ───────────────
fig1, ax1 = plt.subplots(num=1)
ax1.scatter(bx, by, c=colors, cmap="rainbow", s=2, rasterized=True)  # s=2 + rasterized → faster render
ax1.plot(mars_circle.real, mars_circle.imag, "red")
ax1.set(title="B-Plane", xlabel="B_X [km]", ylabel="B_Y [km]")

fig2, ax2 = plt.subplots(num=2)
ax2.scatter(time_, vinf, s=2, rasterized=True)
ax2.set(title="V infinity", xlabel="time [s]", ylabel="V_infinity [km/s]")

fig3, ax3 = plt.subplots(num=3)
ax3.scatter(time_, rsoi, s=2, rasterized=True)
ax3.set(xlabel="time [s]", ylabel="r_soi [km]")

fig4, ax4 = plt.subplots(num=4)
ax4.scatter(time_, bmag, s=2, rasterized=True)  
ax4.set(title="Impact parameter magnitude", xlabel="time [s]", ylabel="Impact Parameter [km]")

if bs_list:
    bs = pd.DataFrame(bs_list, columns=["b_impact_parameter_x", "b_impact_parameter_y"])
    fig5, ax5 = plt.subplots(num=5)
    ax5.scatter(bs["b_impact_parameter_x"], bs["b_impact_parameter_y"])
    ax5.plot(mars_circle.real, mars_circle.imag, "red")
    ax5.set(title="Impact Parameters on B-Plane", xlabel="B_X [km]", ylabel="B_Y [km]")
    ax5.set_aspect('equal')


# ── Animation (on last_df only) ───────────────────────────
# if last_df is not None:
#     df = last_df
#     n  = len(df)

#     # Pre-extract numpy arrays ONCE (iloc in update() is slow per-frame)
#     anim_bx = df["b_impact_parameter_x"].to_numpy()
#     anim_by = df["b_impact_parameter_y"].to_numpy()

#     fig6, ax6 = plt.subplots(figsize=(8, 6))
#     soi_circle = np.exp(1j * np.linspace(0, 2 * np.pi, 100)) * SOI_R
#     ax6.plot(soi_circle.real, soi_circle.imag, "red", label="Target SOI")

#     line,  = ax6.plot([], [], "b-", alpha=0.5)
#     scat   = ax6.scatter([], [], c=[], cmap="rainbow", vmin=0, vmax=n)

#     ax6.set_xlim(anim_bx.min() * 1.1, anim_bx.max() * 1.1)
#     ax6.set_ylim(anim_by.min() * 1.1, anim_by.max() * 1.1)
#     ax6.set(title="B-Plane Trajectory Evolution", xlabel="B_X [km]", ylabel="B_Y [km]")
#     ax6.grid(True, linestyle="--", alpha=0.6)

#     def init():
#         line.set_data([], [])
#         return line,

#     def update(frame):
#         line.set_data(anim_bx[:frame], anim_by[:frame])
#         scat.set_offsets(np.column_stack([anim_bx[:frame], anim_by[:frame]]))
#         scat.set_array(np.arange(frame))
#         return line, scat

#     ani = FuncAnimation(fig6, update, frames=n, init_func=init,
#                         blit=True, interval=20, cache_frame_data=False)  # ← don't cache frames

plt.show()