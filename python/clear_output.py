import shutil
import os

def clear_output():
    directory = 'output/trials'

    # Remove directory and all contents
    if os.path.exists(directory):
        shutil.rmtree(directory)

    # Recreate empty directory
    os.makedirs(directory)

    print(f"Cleared output folder")

clear_output()