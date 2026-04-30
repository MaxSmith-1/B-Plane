import re

input_file = "python/veo_input.txt"   # paste your raw data here
output_file = "python/veo_ref_traj.csv"

rows = []
with open(input_file, "r") as f:
    lines = f.readlines()

# Extract only numeric data lines (skip headers/$$SOE/etc.)
for line in lines:
    if re.match(r'^\s*\d+\.\d+,', line):
        parts = [x.strip() for x in line.split(",")]

        JDTDB = float(parts[0])
        X = float(parts[2])
        Y = float(parts[3])
        Z = float(parts[4])
        VX = float(parts[5])
        VY = float(parts[6])
        VZ = float(parts[7])

        rows.append([JDTDB, X, Y, Z, VX, VY, VZ])

# Convert to time in seconds from first entry
t0 = rows[0][0]
output_rows = []

for i, row in enumerate(rows):
    time_sec = (row[0] - t0) * 86400  # days → seconds
    output_rows.append([int(time_sec)] + row)

# Write CSV
with open(output_file, "w") as f:
    f.write("time,JDTDB,X,Y,Z,VX,VY,VZ\n")
    for r in output_rows:
        f.write(",".join(f"{val:.15E}" if isinstance(val, float) else str(val) for val in r) + "\n")

print(f"Saved to {output_file}")