#include <iostream>
#include <json/json.h>

#include <fstream>
#include <string>

#include <Eigen/Dense>
#include <vector>

#include <ctime>
#include <mutex>



class Simulator
{

    public:

    // Initialze Constructor
    Simulator(double tf, Json::Value spacecraft, Json::Value central_body, bool monte_carlo);
    
    // Initialize simulation loop function
    void simulate(int t_n);

    private:

    // Input variables
    double tf;
    Json::Value spacecraft;
    Json::Value central_body;
    bool monte_carlo;
    // int thread_num = 0;

    double abs_tol;
    double rel_tol;
    // int burn_counter;
    int num_burns;

    // Target body variables
    std::string target_body;
    double target_body_mu;

    // Cached planet coordinates
    std::vector<std::vector<double>> planet_matrix;

    // Gravitational Constant
    double G = 6.674e-20;
    double mu;

 
    // Initialize state vectors
    Eigen::VectorXd state;
    Eigen::VectorXd state_std;
    
    //Eigen::VectorXd derived_state;

    // Initialize time column
    // std::vector<double> time;
    // std::mutex data_mutex;

    // Initialize list for state vector storage
    // std::vector<Eigen::VectorXd> states;
    // std::vector<Eigen::VectorXd> derived_states;

    // Function that gets called on every simulation loop
    void ode_function(const Eigen::VectorXd &x, Eigen::VectorXd &dxdt, const double t);

    // Function that calculates derived state values on each simulation loop
    Eigen::VectorXd build_derived_state(Eigen::VectorXd state, double t);

    // Function that writes states to output csv
    void write_output(std::vector<double>& time, std::vector<Eigen::VectorXd>& states, std::vector<Eigen::VectorXd>& derived_states, int thread_num);

    // Function that writes callback from NASA Horizons API
    static size_t write_callback(void* contents, size_t size, size_t nmemb, std::string* userp);

    // Function that converts unix time to specific datestring for NASA Horizons API
    std::string unix_to_string(int unixTimeSeconds);

    // Function to load planet data into ode function
    void load_planet_data();

    // Function to get target coordinates
    Eigen::Vector3d get_target(double t, Eigen::Vector3d r_sc_rel_sun, bool position);
    
    // Function to generate gaussian state for monte carlo run
    Eigen::VectorXd mc_state(Eigen::VectorXd state);

};