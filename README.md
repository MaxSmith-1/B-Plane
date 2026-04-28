Steps for Compiling on WSL:

0. Download mingw-32 / cmake

1. Download conan: ```pip install conan```

2. Use conan to install dependencies: ```conan install .. --output-folder=. --build=missing --profile=../linux_profile```


2. Compile: ```cmake .. -DCMAKE_TOOLCHAIN_FILE=conan_toolchain.cmake -DCMAKE_BUILD_TYPE=Debug``` then ```make -j``` 

3. Cd into root directory and to run ```./executable/Orbit_Simulator```

Demo Run: ```./executable/B_Plane 18000000 config/Spacecraft/MRO.json config/Bodies/sun.json```
          ```./executable/B_Plane 13000000 config/Spacecraft/VEO.json config/Bodies/sun.json```

Pull spacecraft initial state from: https://ssd.jpl.nasa.gov/horizons/app.html#/
- Use Sun (body center) [500@10] as central body



