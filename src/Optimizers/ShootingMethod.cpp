
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <limits>
#include <thread>


#include <json/json.h>


#include <Optimizers/ShootingMethod.h>
#include <Functions.h>
#include <Simulator.h>

ShootingMethod::ShootingMethod(double tf, double t_burn, Json::Value spacecraft, Json::Value central_body) 
    : Optimizer(tf, t_burn, spacecraft, central_body) {}


Eigen::Vector3d ShootingMethod::optimize(){

    double residual = std::numeric_limits<double>::infinity();
    double threshold = spacecraft["b-threshold"].asDouble();

    double bx_target = spacecraft["target_b_plane_coordinates"]["x"].asDouble();
    double by_target = spacecraft["target_b_plane_coordinates"]["y"].asDouble();

    Eigen::Vector3d velocity(0,0,0);

    double count = 0;
    
    // If no t_burn was specified, optimize for burn time by sweeping Jacobean matricies
    // Get burn time resulution down to the hour

    std::cout << "Optimizing burn time" << std::endl;
    double t_burn_final = 0;     
    if (t_burn == 0.0){

        double t_max = tf;

        // Must be at least one week from planned burn
        double t_min = 86400 * 7;
        
        
        double counter = 0;
        
        while (t_max - t_min > 3600){
            std::cout << "Time trial " << counter << std::endl;
            std::cout << t_max - t_min << counter << std::endl;


            double max_singular_value = -1*std::numeric_limits<double>::infinity();
            unsigned int best_i = 0;

            for(unsigned int i=0; i < 10; i++){

                // Change burn time
                t_burn = t_min + i * ((t_max - t_min) / 10);

                // Add 0 burn to plan
                add_burn(velocity);

                // Compute Jacobean at this time
                Eigen::Matrix<double, 2, 3> J = get_Jacobean(velocity);

                // Get max singular value from J
                Eigen::JacobiSVD<Eigen::Matrix<double, 2, 3>> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
                double singular_value = svd.singularValues()(0);

                // Store the max singular value / time bounds
                if (singular_value > max_singular_value){
                    max_singular_value = singular_value;
                    t_burn_final = t_burn; 
                    best_i = i;
                    
                }

            }

            double current_t_max = t_max;
            double current_t_min = t_min; 
            
            // Shrink the time range and repeat
            t_max = best_i < 9  ? current_t_min + (best_i+1) * ((current_t_max - current_t_min) / 10) : current_t_min + best_i * ((current_t_max - current_t_min) / 10);
            t_min = best_i > 0 ? current_t_min + (best_i-1) * ((current_t_max - current_t_min) / 10) : current_t_min + best_i * ((current_t_max - current_t_min) / 10);

        }

        t_burn = t_burn_final;

    }

    Eigen::Vector3d dv(0,0,0);  


    while(residual > threshold){

        // Clear burn plan / add guess burn
        std::cout << "Adding new burn: " << std::endl;
        add_burn(velocity);

        // Get b-plane coordinates
        std::cout << "Running nominal sim: " << std::endl;
        std::array<double, 2> b_coords = get_b_coordinates(tf, spacecraft, central_body, false);

        int backtrack_count = 0;
        bool cant_pass = false;
        while(std::isnan(b_coords[0]) || (b_coords[0] == 0 && b_coords[1] == 0)){
            std::cout << "WARNING: B-plane not crossed, backtracking" << std::endl;
            if(backtrack_count++ > 10) { 
                std::cerr << "ERROR: Could not find B-plane crossing" << std::endl;
                cant_pass = true;
                break; 
            }
            velocity = velocity - dv * 0.5;
            dv = dv * 0.5;  // shrink dv so next backtrack is smaller
            add_burn(velocity);
            b_coords = get_b_coordinates(tf, spacecraft, central_body, false);
            
        }

        if(cant_pass){
            Eigen::Vector3d return_val(0,0,0);
            return return_val;
        }

        double bx_guess = b_coords[0];
        double by_guess = b_coords[1];

        std::cout << bx_guess << std::endl;

        std::cout << by_guess << std::endl;
        
        // Compute residual vector
        Eigen::Vector2d F(bx_guess - bx_target, by_guess - by_target);

        residual = F.norm();

        std::cout << "Residuals: " << std::endl;
        std::cout << F[0] << std::endl;

        std::cout << F[1] << std::endl;

        // Compute Jacobean (for delta v at t_burn)
        std::cout << "Building Jacobean" << std::endl;
        Eigen::Matrix<double, 2, 3> J = get_Jacobean(velocity);

        std::cout << "Jacobean: " << std::endl;
        std::cout << J << std::endl;


        // Compute next guess / iterate
        Eigen::Matrix<double, 3, 3> JtJ = J.transpose() * J;
        std::cout << "JtJ determinant: " << JtJ.determinant() << std::endl;
        dv = -J.bdcSvd(Eigen::ComputeFullU | Eigen::ComputeFullV).solve(F);
        
        // Velocity is in ICRF coordinate frame
        // TODO: rotate to local RTN frame

        std::cout << "dv: " << std::endl;
        std::cout << dv << std::endl;
        velocity = velocity + dv;

        count += 1;
        
        std::cout << "Trial: " << count << std::endl;
        std::cout << "New residual magnitude:" << std::endl;
        std::cout << F.norm() << std::endl;

        std::cout << "Goal residual" << std::endl;
        std::cout << threshold << std::endl;

        std::cout << "New suggested burn: " << std::endl;
        std::cout << velocity << std::endl;
    }
    

    // Run monte carlo trial with new velocity when done solving

    add_burn(velocity);
    Simulator sim_mc(tf, spacecraft, central_body, true);

    std::cout << "mcing" << std::endl;
    for(int i = 0; i < 64; i+=10){

        //number_trials=64 : number threads=10

        std::vector<std::thread> threads;

        int end = std::min(i + 10, 64);

        for(int j = i; j < end; j++){
        
            threads.push_back(std::thread(&Simulator::simulate, &sim_mc, j));
    
        }

        for (auto& t : threads) {

            t.join();
        }
    
    }

    std::cout << "New suggested burn (ICRF XYZ): " << std::endl;
    std::cout << velocity << std::endl;

    std::cout << "Optimal RTN burn for desired b-plane coordinates at time " << t_burn << " s from start:" << std::endl;

    return get_rtn_burn(velocity, tf, spacecraft, central_body, false);

}


std::array<double, 2> ShootingMethod::get_b_coordinates(double tf, Json::Value spacecraft, Json::Value central_body, bool monte_carlo){

    // Run sim
    sim.simulate(0);

    // Pull impact_parameter_x and impact_parameter_y from output csv
    
    std::string spacecraft_name = spacecraft["name"].asString(); 
    std::ifstream file("output/trials/" + spacecraft_name + "_0.csv");
    
    // Get impact_parameter_x, impact_parameter_y, passed_b_plane columns

    std::map<std::string, int> headers;
    std::string line;
    std::getline(file, line);

    std::stringstream ss(line);
    std::string col;

    int i = 0;

    while (std::getline(ss, col, ',')) {
        headers[col] = i++;
    }

    int bx_col = headers["b_impact_parameter_x"];
    int by_col = headers["b_impact_parameter_y"];
    int passed_b_col = headers["passed_b_plane"];

    double bx_guess = 0;
    double by_guess = 0;
    bool passed_b = false;

    double bx_temp = 0;
    double by_temp = 0;
    bool passed_b_temp = 0; 

    while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string val;
    int col = 0;

        while (std::getline(ss, val, ',')) {
            if (col == bx_col){
                bx_temp = std::stod(val);
            }
            if (col == by_col){
                by_temp = std::stod(val);
            }
            if (col == passed_b_col){ 
                passed_b_temp = std::stoi(val);
            }
            col++;
        }
        //std::cout << passed_b_temp << std::endl;
        if(passed_b_temp == 1){
            bx_guess = bx_temp;
            by_guess = by_temp;
            passed_b = true;
            break;
        }
    }


    if(passed_b){

        return {bx_guess, by_guess};
    }

    else{
            return {0,0};
    }

   
    }

Eigen::Vector3d ShootingMethod::get_rtn_burn(Eigen::Vector3d velocity, double tf, Json::Value spacecraft, Json::Value central_body, bool monte_carlo){


    // Run sim
    sim.simulate(0);

    // Pull impact_parameter_x and impact_parameter_y from output csv
    
    std::string spacecraft_name = spacecraft["name"].asString(); 
    std::ifstream file("output/trials/" + spacecraft_name + "_0.csv");
    
    // Get impact_parameter_x, impact_parameter_y, passed_b_plane columns

    std::map<std::string, int> headers;
    std::string line;
    std::getline(file, line);

    std::stringstream ss(line);
    std::string col;

    int i = 0;

    while (std::getline(ss, col, ',')) {
        headers[col] = i++;
    }

    int burn_col = headers["local_burn_counter"]; 
    Eigen::VectorXd burn_state(6);

    while (std::getline(file, line)) {
    std::stringstream ss(line);
    std::string val;
    int col = 0;

        while (std::getline(ss, val, ',')) {

            if (col < 6){
                burn_state[col] = std::stod(val);
            }
            if (col == burn_col){
                if(std::stod(val) == 1){
                    break;
                }
            }
            col++;
        }
    }


    return Functions::icrf_to_rtn_delta_v(burn_state, velocity);
}



Eigen::Matrix<double, 2, 3> ShootingMethod::get_Jacobean(Eigen::Vector3d velocity){

    double eps = spacecraft["eps"].asDouble();

    Eigen::Matrix<double, 2, 3> J;


    for(unsigned int i = 0; i < 3; i++){

        velocity[i] += eps;
        add_burn(velocity);
        std::array<double, 2> b_coords_plus = get_b_coordinates(tf, spacecraft, central_body, false);

        velocity[i] -= 2*eps;
        add_burn(velocity);
        std::array<double, 2> b_coords_minus = get_b_coordinates(tf, spacecraft, central_body, false);

        double dbx = (b_coords_plus[0] - b_coords_minus[0]) / (2*eps);
        
        double dby = (b_coords_plus[1] - b_coords_minus[1]) / (2*eps);

        J(0, i) = dbx;
        J(1, i) = dby;

        // First good step yay!


        std::cout << "b_coords_plus: " << b_coords_plus[0] << ", " << b_coords_plus[1] << std::endl;
        std::cout << "b_coords_minus: " << b_coords_minus[0] << ", " << b_coords_minus[1] << std::endl;
        std::cout << "dbx: " << dbx << " dby: " << dby << std::endl;
        
        // Reset velocity
        velocity[i] += eps;

    }

    return J;


}

void ShootingMethod::add_burn(Eigen::Vector3d velocity){
    
    // Remove all planned burns
    if (spacecraft.isMember("burns")) {
        spacecraft.removeMember("burns");
    }   


    // Create a single burn entry (the dictionary inside the array)
    Json::Value burn;
    burn["time"] = t_burn;


    std::cout << "t_burn: " << t_burn << std::endl;
    // Create the delta_v_icrf array
    Json::Value delta_v(Json::arrayValue);
    delta_v.append(velocity[0]);
    delta_v.append(velocity[1]);
    delta_v.append(velocity[2]);
    burn["delta_v_icrf"] = delta_v;

    // Create the burns array and append the burn entry
    Json::Value burns(Json::arrayValue);
    burns.append(burn);

    // Add to your parent object
    spacecraft["burns"] = burns;


    sim.set_spacecraft(spacecraft);

    // std::cout << "Burns in spacecraft after add_burn: " << spacecraft["burns"] << std::endl;
    // std::cout << "Delta v set to: " << velocity[0] << ", " << velocity[1] << ", " << velocity[2] << std::endl;


    
}

