#include <iostream>

#include <json/json.h>
#include <fstream>
#include <string>
#include <thread>

#include <Simulator.h>


int main(int argc, char* argv[]){


    std::cout << "Welcome to the B-Plane Optimizer" <<std::endl;

    // Read in inputs
    int tf = std::stoi(argv[1]);
    std::string spacecraft = argv[2];
    std::string central_body = argv[3];

    int mc_trials = 0;
    int num_threads = 1;

    if(argc > 4){

        mc_trials = std::stoi(argv[4]);
        num_threads = std::stoi(argv[5]);
        std::cout << "Monte Carlo run" << std::endl;

    }

    // argv[6,7] = Target b-plane x,y coordinates
    // TODO: Optimizer class that calls simulator 


    
    Json::Value spacecraft_json, body_json;
    
    // Load data from json files
    try{
        std::ifstream input1(spacecraft);
        std::ifstream input2(central_body);

        input1 >> spacecraft_json;
        input2 >> body_json;
    }
    catch(const std::exception& e){

        std::cout << "Could not load json files" << std::endl;

        std::cout << e.what() << std::endl;
        return 1;
    }
    
    // Call simulator class

    bool mc = mc_trials > 0;
    std::cout << "Simulating spacecraft " << spacecraft_json["name"].asString() << " for " << std::to_string(tf) << "s." << std::endl;
    Simulator sim(tf, spacecraft_json, body_json, mc);


    // Monte Carlo run with multithreading
    if(mc){

        for(int i = 0; i < mc_trials; i+=num_threads){

            std::vector<std::thread> threads;

            int end = std::min(i + num_threads, mc_trials);

            for(int j = i; j < end; j++){
            
                threads.push_back(std::thread(&Simulator::simulate, &sim, j));
        
            }

            for (auto& t : threads) {

                t.join();
            }
        
        }

    }

    // One off run
    else{
        sim.simulate(0);
        
    }
    // Call this function in python and generate some basic plots

    return 0;
}