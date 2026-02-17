#include "Simulator.h"

#include <iostream>
#include <json/json.h>
#include <fstream>
#include <string>
#include <cmath>
#include <random>

#include <filesystem>
#include <sys/stat.h>

#include <iomanip>

//#define BOOST_MATH_DISABLE_THREADS  // Disable threading in Boost.Math

#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/stepper/generation/generation_runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>
#include <Functions.h>

#include <ctime>
#include <curl/curl.h>

#include <sstream>
#include <ctime>
#include <chrono>
#include <typeinfo>


using namespace boost::numeric::odeint;


Simulator::Simulator(double tf, Json::Value spacecraft, Json::Value central_body, bool monte_carlo){


    this->tf = tf;
    this->spacecraft = spacecraft;
    this->central_body = central_body;
    this-> monte_carlo = monte_carlo;


    mu = G*central_body["mass"].asDouble();

    // Set up initial state
    state.resize(6);

    state << spacecraft["initial_condition"]["ICRF_position"][0].asDouble(),
             spacecraft["initial_condition"]["ICRF_position"][1].asDouble(),
             spacecraft["initial_condition"]["ICRF_position"][2].asDouble(),
             spacecraft["initial_condition"]["ICRF_velocity"][0].asDouble(),
             spacecraft["initial_condition"]["ICRF_velocity"][1].asDouble(),
             spacecraft["initial_condition"]["ICRF_velocity"][2].asDouble();


    // If a monte carlo run, perturb the state
    abs_tol = spacecraft["abs_tol"].asDouble();
    rel_tol = spacecraft["rel_tol"].asDouble();

    num_burns = spacecraft["burns"].size();

    // TODO: Set this up for spacecraft["target_body"]

    // Set up settings

    Json::Value target = spacecraft["target_body"];
    Json::Value target_mass = spacecraft["target_body_mass"];

    target_body = target.asString();
    target_body_mu = G * target_mass.asDouble();



    // Planetary tol is the number of divisions
    std::string planetary_step_size = std::to_string(static_cast<int>(tf * spacecraft["planet_tol"].asDouble()));

    // Dates are inputted as unix time, outputted as date strings
    std::string initial_date = unix_to_string(spacecraft["date"].asInt());
    std::string final_date = unix_to_string(spacecraft["date"].asInt() + tf);

    
    
    // Initialize libcurl
    CURL* curl;
    CURLcode res;
    std::string responseData;
    curl_global_init(CURL_GLOBAL_DEFAULT);
    curl = curl_easy_init();


    if(curl) {

        std::string input_url = "https://ssd.jpl.nasa.gov/api/horizons.api?format=text";
        input_url += "&COMMAND='" + target_body + "'";
        input_url += "&OBJ_DATA='YES'";
        input_url += "&MAKE_EPHEM='YES'";
        input_url += "&EPHEM_TYPE='VECTORS'";
        input_url += "&CSV_FORMAT='YES'";
        input_url += "&CENTER='500@10'";
        input_url += "&VEC_TABLE='2'"; // Position and Velocity state
        input_url += "&START_TIME='" + initial_date + "'";
        input_url += "&STOP_TIME='" + final_date + "'";
        input_url += "&STEP_SIZE='" + planetary_step_size + "'";

        std::cout << input_url;
        curl_easy_setopt(curl, CURLOPT_URL, input_url.c_str());
        
        // Set the callback function to handle response
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
        curl_easy_setopt(curl, CURLOPT_WRITEDATA, &responseData);
        
        // Follow redirects
        curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
        
        // Perform the request
        res = curl_easy_perform(curl);

        //std::cout << res << std::endl;
        
        // Check for errors
        if(res != CURLE_OK) {
            std::cerr << "curl_easy_perform() failed: " 
                    << curl_easy_strerror(res) << std::endl;
        } else {
            std::cout << "Response: " << responseData << std::endl;
        }
        
        //std::cout << typeid(responseData).name() << std::endl;
        
        // Cleanup
        curl_easy_cleanup(curl);

        // Parse the data for xyz coordinates
        size_t startPos = responseData.find("$$SOE");

        
        startPos += 5;  // Move past the start delimiter
        
        size_t endPos = responseData.find("$$EOE", startPos);
        
        std::string body_data = responseData.substr(startPos, endPos - startPos);

        // std::cout << body_data << std::endl;

        std::string output_string = "output/planets/" + target_body + ".csv";

        std::ofstream file(output_string);

        file << body_data;
        file.close();

    }
    curl_global_cleanup();

    // Load in planet coordinates
    load_planet_data();


    }
    

void Simulator::simulate(int t_n){

    // time / state / derived state variables changed to be local to this function to prevent overwriting when multithreading
    std::vector<double> local_time;
    std::vector<Eigen::VectorXd> local_states;
    std::vector<Eigen::VectorXd> local_derived_states;
    Eigen::VectorXd local_state = state;  // Copy initial state

    if(monte_carlo){
        std::cout << "Monte Carloing" << std::endl;
        state = mc_state(state);
        std::cout << state[0] << std::endl;
    }


    int local_burn_counter = 0;

    // Define rk45 solver
    typedef runge_kutta_dopri5<Eigen::VectorXd, double, Eigen::VectorXd, double, vector_space_algebra> error_stepper_type;

    auto ode_func = [this](const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double t) {
        this->ode_function(x, dxdt, t);
    };

    // Observer lambda function that builds derived state / checks for burns every time step
    auto observer = [&](Eigen::VectorXd &state_ref, double t) {
        Eigen::VectorXd derived_state = build_derived_state(state_ref, t);
        
        local_derived_states.push_back(derived_state);
        local_states.push_back(state_ref);
        local_time.push_back(t);

        // Check for collision
        if (std::hypot(state_ref[0], state_ref[1], state_ref[2]) < central_body["radius"].asDouble()){
            throw std::runtime_error("Spacecraft impacted the central body surface.");
        }

        // Check for burns
        if(t >= spacecraft["burns"][local_burn_counter]["time"].asDouble() && local_burn_counter < num_burns){
            std::cout << "Executing burn" << std::endl;
            state_ref[3] += spacecraft["burns"][local_burn_counter]["delta_v_icrf"][0].asDouble();
            state_ref[4] += spacecraft["burns"][local_burn_counter]["delta_v_icrf"][1].asDouble();
            state_ref[5] += spacecraft["burns"][local_burn_counter]["delta_v_icrf"][2].asDouble();
            local_burn_counter++;
        }
    };

    try{
        integrate_adaptive(make_controlled<error_stepper_type>(abs_tol, rel_tol), 
                          ode_func, local_state, 0.0, tf, 0.01, observer);
    }
    catch(const std::runtime_error& e){
        std::cout << "Spacecraft impacted the central body surface." << std::endl;
    }

    // Write outputs
    write_output(local_time, local_states, local_derived_states, t_n);
}

// Function that gets called on every simulation loop
void Simulator::ode_function(const Eigen::VectorXd &x, Eigen::VectorXd &dxdt, const double t){

    // Initialize derivitive vector
    if (dxdt.size() != x.size()) {
        dxdt.resize(x.size());
    }
    dxdt.setZero();

    Eigen::Vector3d r_sc_rel_sun = x.head(3);

    // Add gravitational effect of target body (just need the coordinates)
    // ASSUME just 1 target body
    Eigen::Vector3d r_sc_rel_target = get_target(t, r_sc_rel_sun, true);

    // std::cout << r_target_rel_sun[0] << std::endl;
    
    double target_accel_x = (-target_body_mu / pow(r_sc_rel_target.norm(), 3)) * r_sc_rel_target[0];
    double target_accel_y = (-target_body_mu / pow(r_sc_rel_target.norm(), 3)) * r_sc_rel_target[1];
    double target_accel_z = (-target_body_mu / pow(r_sc_rel_target.norm(), 3)) * r_sc_rel_target[2];

    // Two Body Problem + Gravity of Target Body
    dxdt[0] = x[3];
    dxdt[1] = x[4];
    dxdt[2] = x[5];

    dxdt[3] = (-(mu / std::pow(r_sc_rel_sun.norm(), 3)) * r_sc_rel_sun[0]) + target_accel_x;
    dxdt[4] = (-(mu / std::pow(r_sc_rel_sun.norm(), 3)) * r_sc_rel_sun[1]) + target_accel_y;
    dxdt[5] = (-(mu / std::pow(r_sc_rel_sun.norm(), 3)) * r_sc_rel_sun[2]) + target_accel_z;

}

// Function that calculates derived state values on each simulation loop
Eigen::VectorXd Simulator::build_derived_state(Eigen::VectorXd state, double t){

    double rad_to_deg = 180 / M_PI;
    // DERIVED STATES TO CALCULATE
    // a,e,i,raan, omega, f, E, M, n, p, h, flight path angle

    // B-Plane variables
    // r_rel_target, v_rel_target |b|, S, T, R, B vectors
    Eigen::Vector3d r_vec(state[0], state[1], state[2]);
    Eigen::Vector3d v_vec(state[3], state[4], state[5]);

    double r = r_vec.norm();  
    double v = v_vec.norm();
    // Specific energy (vis-viva) equation 
    double Energy = (std::pow(v, 2) / 2) - (mu / r);

    // Semi-major axis
    double a = -mu / (2*Energy);
    
    // Mean motion
    double n = std::sqrt(mu / pow(a, 3));

    // Orbital Period
    double T = 2*M_PI/n;

    // Angular Momentum Vector
    Eigen::Vector3d h_vec = r_vec.cross(v_vec);
    double h = h_vec.norm();

    // Eccentricity Vector
    Eigen::Vector3d e_vec = ((1/mu) * (v_vec.cross(h_vec))) - (r_vec/r);
    double e = e_vec.norm();

    // Semi-latus rectrum
    double p = std::pow(h, 2) / mu;

    // Semi-minor axis
    double b = a*std::sqrt(1 - pow(e, 2));

    // Radii of apogee and perigee
    double ra = a*(1+e);
    double rp = a*(1-e);

    // True Anomaly
    double f = std::acos(((p / r) - 1) / e);

    if (r_vec.dot(v_vec) < 0) {
        f = 2.0 * M_PI - f;
    }
    // Eccentric Anomaly
    double E = 2*std::atan(std::sqrt((1-e) / (1+e))*std::tan(f/2));

    // Mean Anomaly
    double M = E - e*std::sin(E);

    // Flight Path Angle
    double gamma = std::atan((e*std::sin(f)) / (1+ e*std::cos(f)));

    // Inclination
    Eigen::Vector3d x(1, 0, 0);
    Eigen::Vector3d y(0, 1, 0);
    Eigen::Vector3d z(0, 0, 1);

    double i = std::acos((h_vec / h).dot(z));

    // Longitude of Ascending Node
    double laan = std::acos(x.dot(z.cross(h_vec / h)) / std::sin(i));

    if(laan < 0){
        laan += 2*M_PI;
    }

    // Argument of Periapsis
    Eigen::Vector3d node_vector = z.cross(h_vec);
    double omega = std::acos(node_vector.dot(e_vec) / (e*node_vector.norm())); 
    
    if (e_vec.z() < 0){
        omega = 2.0 * M_PI - omega;
    }



    // Should probably spit out sphere of influence



    // Build B-Plane variables
    // r_rel_target, v_rel_target |b|, S, T, R, B vectors
    Eigen::Vector3d r_sc_rel_target = get_target(t, r_vec, true);
    Eigen::Vector3d v_sc_rel_target = get_target(t, v_vec, false);
    Eigen::Vector3d h_sc_rel_target = r_sc_rel_target.cross(v_sc_rel_target);
    Eigen::Vector3d e_sc_rel_target = ((1/target_body_mu) * (v_sc_rel_target.cross(h_sc_rel_target))) - (r_sc_rel_target/r_sc_rel_target.norm());

    double r_soi = std::pow((1.0 * target_body_mu / mu), 0.4) * a;

    bool in_target_soi = r_sc_rel_target.norm() < r_soi;

    double V_infinity = std::sqrt(std::pow(v_sc_rel_target.norm(), 2) - ((2*target_body_mu) / r_sc_rel_target.norm()));
    double b_impact_parameter = (h_sc_rel_target).norm() / V_infinity;
    double beta = std::acos(1.0 / e_sc_rel_target.norm());

    // B-Plane coordinate frame vectors     
    Eigen::Vector3d S_hat = (std::cos(beta) / e_sc_rel_target.norm())*e_sc_rel_target + (std::sin(beta) / (e_sc_rel_target.norm() * h_sc_rel_target.norm())) * h_sc_rel_target.cross(e_sc_rel_target);
    Eigen::Vector3d T_hat = S_hat.cross(Eigen::Vector3d::UnitZ()).normalized();
    Eigen::Vector3d R_hat = S_hat.cross(T_hat);

    // Impact parameter vector
    Eigen::Vector3d B_vec = b_impact_parameter * (S_hat.cross(h_sc_rel_target).normalized());

    double b_impact_parameter_x =  B_vec.dot(T_hat);
    double b_impact_parameter_y = B_vec.dot(R_hat);

    // When this quantity goes positive, spacecraft has passed through B-plane, take b_x, y at this point to be impact parameter
    double b_impact_parameter_s_component = r_sc_rel_target.dot(S_hat);

    bool passed_b_plane = b_impact_parameter_s_component > 0;

    


    // Convert angular quantities to degrees
    f = f * rad_to_deg; 
    E = E * rad_to_deg; 
    M = M * rad_to_deg; 
    gamma = gamma * rad_to_deg; 
    i = i * rad_to_deg; 
    laan = laan * rad_to_deg; 
    omega = omega * rad_to_deg; 


    Eigen::VectorXd derived_state(34);

    derived_state << v, r, Energy, a, n, T, h, h_vec[0], h_vec[1], h_vec[2], 
                 e, e_vec[0], e_vec[1], e_vec[2], p, ra, rp, 
                 b, f, E, M, gamma, i, laan, omega,
                 V_infinity, b_impact_parameter, beta, b_impact_parameter_x, b_impact_parameter_y,
                 r_soi, in_target_soi,
                 b_impact_parameter_s_component, passed_b_plane;
 
    return derived_state;

}

// Function that writes states to output csv
void Simulator::write_output(std::vector<double>& time, 
                             std::vector<Eigen::VectorXd>& states, 
                             std::vector<Eigen::VectorXd>& derived_states,
                             int thread_num){
    
    // Filename
    const std::string& filename = spacecraft["name"].asString() + "_" + std::to_string(thread_num) + ".csv";

    // Hardcoded headers
    const std::vector<std::string> state_headers = {
        "ICRF_X", "ICRF_Y", "ICRF_Z", "Vx", "Vy", "Vz"
    };
    
    const std::vector<std::string> derived_headers = {
        "v", "r", "E", "a", "n", "T", "h", "h_x", "h_y", "h_z", 
        "e", "e_x", "e_y", "e_z", "p", "ra", "rp", 
        "b", "f", "E_anom", "M", "gamma", "i", "laan", "omega",
        "V_infinity", "b_impact_parameter", "beta", "b_impact_parameter_x", "b_impact_parameter_y",
        "r_soi", "within_target_soi", "b_impact_parameter_s_component", "passed_b_plane"
    };
    
    // Create output directory if it doesn't exist
    std::string output_dir = "output/trials";

    std::filesystem::create_directories(output_dir.c_str());


    // Construct full file path,
    std::string filepath = output_dir + "/" + filename;
    
    std::ofstream file(filepath);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filepath << std::endl;
        return;
    }
    
    // Set precision
    file << std::scientific << std::setprecision(10);
    
    // Check that all vectors have the same length

    std::cout << time.size() << std::endl;
    std::cout << states.size() << std::endl;
    std::cout << derived_states.size() << std::endl;


    if (time.size() != states.size() || time.size() != derived_states.size()) {
        std::cerr << "Error: time, states, and derived_states must have the same length" << std::endl;
        file.close();
        return;
    }
    
    if (time.empty()) {
        std::cerr << "Warning: No data to write" << std::endl;
        file.close();
        return;
    }
    
    // Get dimensions
    int num_states = states[0].size();
    int num_derived = derived_states[0].size();
    
    // Write header row
    file << "time";
    
    // State headers
    for (unsigned int i = 0; i < num_states; i++) {
        file << ",";
        if (i < state_headers.size()) {
            file << state_headers[i];
        } else {
            file << "state_" << i;
        }
    }
    
    // Derived state headers
    for (unsigned int i = 0; i < num_derived; i++) {
        file << ",";
        if (i < derived_headers.size()) {
            file << derived_headers[i];
        } else {
            file << "derived_" << i;
        }
    }
    file << "\n";
    
    // Write data rows
    for (size_t row = 0; row < time.size(); row++) {
        // Write time
        file << time[row];
        
        // Write state values
        for (int i = 0; i < states[row].size(); i++) {
            file << "," << states[row](i);
        }
        
        // Write derived state values
        for (int i = 0; i < derived_states[row].size(); i++) {
            file << "," << derived_states[row](i);
        }
        
        file << "\n";
    }
    
    file.close();
    std::cout << "Data written to " << filepath << std::endl;
}


size_t Simulator::write_callback(void* contents, size_t size, size_t nmemb, std::string* userp) {
    size_t totalSize = size * nmemb;
    userp->append((char*)contents, totalSize);
    return totalSize;
}

std::string Simulator::unix_to_string(int unixTimeSeconds) {
    // Convert to time_t
    std::time_t time = static_cast<std::time_t>(unixTimeSeconds);
    
    // Convert to tm structure (UTC)
    std::tm* tm = std::gmtime(&time);
    
    // Month abbreviations
    const char* months[] = {
        "01", "02", "03", "04", "05", "06",
        "07", "08", "09", "10", "11", "12"
    };
    
    // Format the string
    std::ostringstream oss;
    oss << std::setfill('0')
        << std::setw(4) << (tm->tm_year + 1900) << "-"
        << months[tm->tm_mon] << "-"
        << std::setw(2) << tm->tm_mday << "%20"
        << std::setw(2) << tm->tm_hour << ":"
        << std::setw(2) << tm->tm_min << ":"
        << std::setw(2) << tm->tm_sec << "."
        << std::setw(3) << 0;

    return oss.str();
}

void Simulator::load_planet_data(){
    std::string filename = "output/planets/" + target_body + ".csv";
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "ERROR: Could not open " << filename << std::endl;
        return;
    }
    
    std::string line;
    
    while(std::getline(file, line)){
        std::stringstream ss(line);
        std::string cell_str;
        std::vector<double> row;
        
        while (std::getline(ss, cell_str, ',')) {
            try {
                row.push_back(std::stod(cell_str));
            } catch (const std::invalid_argument& e) {
                // Skip non-numeric cells (like the date string)
                continue;
            }
        }
        
        // Only add rows that have the expected number of values (4: JD, X, Y, Z)
        if (row.size() == 7) {
            planet_matrix.push_back(row);
        }
    }
    file.close();
    
    std::cout << "Loaded " << planet_matrix.size() << " planet coordinate rows" << std::endl;
}

Eigen::Vector3d Simulator::get_target(double t, Eigen::Vector3d sc_rel_sun, bool position){

    // Get x, y, z velocity
    if (position){
        Eigen::Vector3d r_target_rel_sun(0, 0, 0);
        double time_percent = t / tf;

        double planet_t0 = planet_matrix[1][0];

        double planet_tf = planet_matrix[planet_matrix.size() - 1][0];

        double planet_t = planet_t0 + ((planet_tf - planet_t0) * time_percent);

        auto it = std::lower_bound(planet_matrix.begin(), planet_matrix.end(), planet_t,
            [](const std::vector<double>& row, double val) {
                return row[0] < val;
            });

        int n = static_cast<int>(planet_matrix.size());
        int j = std::distance(planet_matrix.begin(), it);

        j = std::max(1, std::min(j, n - 1));

        double mutiplier = (planet_t - planet_matrix[j-1][0]) / (planet_matrix[j][0] - planet_matrix[j-1][0]);

        r_target_rel_sun[0] = (mutiplier * (planet_matrix[j][1] - planet_matrix[j - 1][1])) + planet_matrix[j - 1][1];
        r_target_rel_sun[1] = (mutiplier * (planet_matrix[j][2] - planet_matrix[j - 1][2])) + planet_matrix[j - 1][2];
        r_target_rel_sun[2] = (mutiplier * (planet_matrix[j][3] - planet_matrix[j - 1][3])) + planet_matrix[j - 1][3];


        return sc_rel_sun - r_target_rel_sun;
    }

    // Get x, y, z velocity
    else{
        Eigen::Vector3d v_target_rel_sun(0, 0, 0);
        double time_percent = t / tf;

        double planet_t0 = planet_matrix[1][0];

        double planet_tf = planet_matrix[planet_matrix.size() - 1][0];

        double planet_t = planet_t0 + ((planet_tf - planet_t0) * time_percent);

        auto it = std::lower_bound(planet_matrix.begin(), planet_matrix.end(), planet_t,
            [](const std::vector<double>& row, double val) {
                return row[0] < val;
            });

        int n = static_cast<int>(planet_matrix.size());
        int j = std::distance(planet_matrix.begin(), it);

        j = std::max(1, std::min(j, n - 1));

        double mutiplier = (planet_t - planet_matrix[j-1][0]) / (planet_matrix[j][0] - planet_matrix[j-1][0]);

        v_target_rel_sun[0] = (mutiplier * (planet_matrix[j][4] - planet_matrix[j - 1][4])) + planet_matrix[j - 1][4];
        v_target_rel_sun[1] = (mutiplier * (planet_matrix[j][5] - planet_matrix[j - 1][5])) + planet_matrix[j - 1][5];
        v_target_rel_sun[2] = (mutiplier * (planet_matrix[j][6] - planet_matrix[j - 1][6])) + planet_matrix[j - 1][6];


        return sc_rel_sun - v_target_rel_sun;
    }
  
}

Eigen::VectorXd Simulator::mc_state(Eigen::VectorXd state){


    // Initialize random class
    std::random_device rd;
    std::mt19937 gen(rd());

    state_std.resize(6);
    state_std << spacecraft["std"]["ICRF_position"][0].asDouble(),
            spacecraft["std"]["ICRF_position"][1].asDouble(),
            spacecraft["std"]["ICRF_position"][2].asDouble(),
            spacecraft["std"]["ICRF_velocity"][0].asDouble(),
            spacecraft["std"]["ICRF_velocity"][1].asDouble(),
            spacecraft["std"]["ICRF_velocity"][2].asDouble();
    
    // Perturb initial state into gaussian distribution
    for(unsigned int i = 0; i < state.size(); i++){

        std::normal_distribution<double> dist(state[i], state_std[i]);
        state[i] = dist(gen);
    }
    return state;

}