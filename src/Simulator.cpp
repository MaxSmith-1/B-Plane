#include "Simulator.h"

#include <iostream>
#include <json/json.h>
#include <fstream>
#include <string>
#include <cmath>

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


Simulator::Simulator(double tf, Json::Value spacecraft, Json::Value central_body){


    this->tf = tf;
    this->spacecraft = spacecraft;
    this->central_body = central_body;

    mu = G*central_body["mass"].asDouble();

    // Set up initial state
    state.resize(6);

    state << spacecraft["initial_condition"]["ICRF_position"][0].asDouble(),
             spacecraft["initial_condition"]["ICRF_position"][1].asDouble(),
             spacecraft["initial_condition"]["ICRF_position"][2].asDouble(),
             spacecraft["initial_condition"]["ICRF_velocity"][0].asDouble(),
             spacecraft["initial_condition"]["ICRF_velocity"][1].asDouble(),
             spacecraft["initial_condition"]["ICRF_velocity"][2].asDouble();

    
    abs_tol = spacecraft["abs_tol"].asDouble();
    rel_tol = spacecraft["rel_tol"].asDouble();

    burn_counter = 0;
    num_burns = spacecraft["burns"].size();


    // Set up settings
    // alternate_bodies = {"10"};
    // alternate_bodies_mu = {1.327e11};

    // Json::Value planets = central_body["alternate_bodies"];
    // Json::Value mus = central_body["alternate_bodies_mu"];

    // std::cout << planets << std::endl;

    // for(unsigned int i = 0; i < planets.size(); i++){
    //     alternate_bodies.push_back(planets[i].asString());
    //     alternate_bodies_mu.push_back(mus[i].asDouble());
        
    // }

    // std::string planetary_step_size = std::to_string(static_cast<int>(tf * spacecraft["planet_tol"].asDouble()));
    // std::string initial_date = unix_to_string(spacecraft["date"].asInt());
    // std::string final_date = unix_to_string(spacecraft["date"].asInt() + tf);
    
    // for(unsigned int i = 0; i < alternate_bodies.size(); i++){

    //     // Initialize libcurl
    //     CURL* curl;
    //     CURLcode res;
    //     std::string responseData;
    //     curl_global_init(CURL_GLOBAL_DEFAULT);
    //     curl = curl_easy_init();


    //     if(curl) {

    //         std::string input_url = "https://ssd.jpl.nasa.gov/api/horizons.api?format=text";
    //         input_url += "&COMMAND='" + alternate_bodies[i] + "'";
    //         input_url += "&OBJ_DATA='YES'";
    //         input_url += "&MAKE_EPHEM='YES'";
    //         input_url += "&EPHEM_TYPE='VECTORS'";
    //         input_url += "&CSV_FORMAT='YES'";
    //         input_url += "&CENTER='500@399'";
    //         input_url += "&VEC_TABLE='1'";
    //         input_url += "&START_TIME='" + initial_date + "'";
    //         input_url += "&STOP_TIME='" + final_date + "'";
    //         input_url += "&STEP_SIZE='" + planetary_step_size + "'";

    //         std::cout << input_url;
    //         curl_easy_setopt(curl, CURLOPT_URL, input_url.c_str());
            
    //         // Set the callback function to handle response
    //         curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);
    //         curl_easy_setopt(curl, CURLOPT_WRITEDATA, &responseData);
            
    //         // Follow redirects
    //         curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
            
    //         // Perform the request
    //         res = curl_easy_perform(curl);

    //         //std::cout << res << std::endl;
            
    //         // Check for errors
    //         if(res != CURLE_OK) {
    //             std::cerr << "curl_easy_perform() failed: " 
    //                     << curl_easy_strerror(res) << std::endl;
    //         } else {
    //             std::cout << "Response: " << responseData << std::endl;
    //         }
            
    //         //std::cout << typeid(responseData).name() << std::endl;
            
    //         // Cleanup
    //         curl_easy_cleanup(curl);

    //         // Parse the data for xyz coordinates
    //         size_t startPos = responseData.find("$$SOE");

            
    //         startPos += 5;  // Move past the start delimiter
            
    //         size_t endPos = responseData.find("$$EOE", startPos);
            
    //         std::string body_data = responseData.substr(startPos, endPos - startPos);

    //         // std::cout << body_data << std::endl;

    //         std::string output_string = "output/planets/" + alternate_bodies[i] + ".csv";

    //         if(i == 0){
    //             output_string = "output/planets/sun_coordinates.csv";
    //         }

    //         std::ofstream file(output_string);

    //         file << body_data;
    //         file.close();

    //     }

    //     curl_global_cleanup();
    // }
        

    // #ifdef _WIN32
    //     this->year_start_unix = static_cast<double>(_mkgmtime(&start_of_year));
    // #else
    //     this->year_start_unix = static_cast<double>(timegm(&start_of_year));
    // #endif

    
}

void Simulator::simulate(){

    // Define rk45 solver to work with Eigen library vectors
    typedef runge_kutta_dopri5<Eigen::VectorXd, double, Eigen::VectorXd, double, vector_space_algebra> error_stepper_type;

    // Wrap the ode function in a lambda
    auto ode_func = [this](const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, const double t) {
        this->ode_function(x, dxdt, t);
    };


    // Integrate the ode function until tf or until collision event (detected in observe)
    // Also record states vs. time through observe function
    try{
        integrate_adaptive( make_controlled< error_stepper_type >( abs_tol , rel_tol ) , 
                    ode_func , state , 0.0 , tf , 0.01, [this](Eigen::VectorXd &state , double t ) {
                                this->observe( state , t );
                            });
                        }
    
    catch(const std::runtime_error& e){

        std::cout << "Spacecraft impacted the central body surface." << std::endl;

    }

    // Wtite outputs to csv
    write_output(time, states, derived_states);
}

void Simulator::observe(Eigen::VectorXd &state , double t){
    
    // Build derived state from current state
    Eigen::VectorXd derived_state = build_derived_state(state, t);

    // Save current time, state, and derived states into respective vectors 
    derived_states.push_back(derived_state);
    states.push_back( state );
    time.push_back( t );

    // Check for collision event
    if (std::hypot(state[0], state[1], state[2]) < central_body["radius"].asDouble()){
        throw std::runtime_error("Spacecraft impacted the central body surface.");
    }

    // Check for burn events
    if(t >= spacecraft["burns"][burn_counter]["time"].asDouble() && burn_counter < num_burns){

        std::cout << "Execuring burn at [" << spacecraft["burns"][burn_counter]["delta_v_rtn"][0].asDouble() <<
        ", " << spacecraft["burns"][burn_counter]["delta_v_rtn"][1].asDouble() << ", " <<
        spacecraft["burns"][burn_counter]["delta_v_rtn"][2].asDouble() << "] (rtn) km/s at t = " << 
        t << std::endl;

        // Convert burn in rtn frame to eci frame
        std::vector<double> delta_v_rtn_vec = {
        spacecraft["burns"][burn_counter]["delta_v_rtn"][0].asDouble(),
        spacecraft["burns"][burn_counter]["delta_v_rtn"][1].asDouble(),
        spacecraft["burns"][burn_counter]["delta_v_rtn"][2].asDouble()
        };

        Eigen::Vector3d delta_v_eci(Functions::rtn_to_eci_delta_v(state, delta_v_rtn_vec));
        
        // Add delta-v to current velocity state
        state[3] += delta_v_eci[0];
        state[4] += delta_v_eci[1];
        state[5] += delta_v_eci[2];

        burn_counter++;
    }
}

// Function that gets called on every simulation loop
void Simulator::ode_function(const Eigen::VectorXd &x, Eigen::VectorXd &dxdt, const double t){

    // Initialize derivitive vector
    dxdt = Eigen::VectorXd(x) * 0.0;

    double r = x.head(3).norm();

    // Two Body Problem
    dxdt[0] = x[3];
    dxdt[1] = x[4];
    dxdt[2] = x[5];

    dxdt[3] = (-(mu / std::pow(r, 3)) * x[0]);
    dxdt[4] = (-(mu / std::pow(r, 3)) * x[1]);
    dxdt[5] = (-(mu / std::pow(r, 3)) * x[2]);

}

// Function that calculates derived state values on each simulation loop
Eigen::VectorXd Simulator::build_derived_state(Eigen::VectorXd state, double t){

    double rad_to_deg = 180 / M_PI;
    // DERIVED STATES TO CALCULATE
    // a,e,i,raan, omega, f, E, M, n, p, h, flight path angle
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


    // Convert angular quantities to degrees
    f = f * rad_to_deg; 
    E = E * rad_to_deg; 
    M = M * rad_to_deg; 
    gamma = gamma * rad_to_deg; 
    i = i * rad_to_deg; 
    laan = laan * rad_to_deg; 
    omega = omega * rad_to_deg; 

    Eigen::VectorXd derived_state(25);

    derived_state << v, r, Energy, a, n, T, h, h_vec[0], h_vec[1], h_vec[2], 
                 e, e_vec[0], e_vec[1], e_vec[2], p, ra, rp, 
                 b, f, E, M, gamma, i, laan, omega;
 
    return derived_state;

}

// Function that writes states to output csv
void Simulator::write_output(std::vector<double>& time, 
                  std::vector<Eigen::VectorXd>& states, 
                  std::vector<Eigen::VectorXd>& derived_states
                  ) {
    
    // Filename
    const std::string& filename = spacecraft["name"].asString() + "_output.csv";

    // Hardcoded headers
    const std::vector<std::string> state_headers = {
        "ICRF_X", "ICRF_Y", "ICRF_Z", "Vx", "Vy", "Vz"
    };
    
    const std::vector<std::string> derived_headers = {
        "v", "r", "E", "a", "n", "T", "h", "h_x", "h_y", "h_z", 
        "e", "e_x", "e_y", "e_z", "p", "ra", "rp", 
        "b", "f", "E_anom", "M", "gamma", "i", "laan", "omega"
    };
    
    // Create output directory if it doesn't exist
    std::string output_dir = "output";

    std::filesystem::create_directories(output_dir.c_str());


    // Construct full file path
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
    for (int i = 0; i < num_states; i++) {
        file << ",";
        if (i < state_headers.size()) {
            file << state_headers[i];
        } else {
            file << "state_" << i;
        }
    }
    
    // Derived state headers
    for (int i = 0; i < num_derived; i++) {
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
