#include <iostream>
#include <json/json.h>

#include <fstream>
#include <string>

#include <Eigen/Dense>
#include <vector>

#include "Nrlmsise00.hpp"
#include <ctime>



class Simulator
{

    public:

    // Initialze Constructor
    Simulator(double tf, Json::Value spacecraft, Json::Value central_body);
    
    // Initialize simulation loop function
    void simulate();

    private:

    // Input variables
    double tf;
    Json::Value spacecraft;
    Json::Value central_body;

    double abs_tol;
    double rel_tol;
    int burn_counter;
    int num_burns;

    // Gravitational Constant
    double G = 6.674e-20;
    double mu;

    // Time Constants
    time_t total_seconds;
    struct tm *now;
    struct tm start_of_year;
    double year_start_unix;

    // Drag class
    std::array<int, 24> atmos_flags; 
    atmos::CNrlmsise00 nrlmsise00;

    // Initialize state vectors
    Eigen::VectorXd state;
    //Eigen::VectorXd derived_state;

    // Initialize time column
    std::vector<double> time;

    // Initialize list for state vector storage
    std::vector<Eigen::VectorXd> states;
    std::vector<Eigen::VectorXd> derived_states;

    // Initialize list of gravitational bodies
    std::vector<std::string> alternate_bodies;
    std::vector<double> alternate_bodies_mu;

    // Function that gets called on every simulation loop
    void ode_function(const Eigen::VectorXd &x, Eigen::VectorXd &dxdt, const double t);

    // Function that calculates derived state values on each simulation loop
    Eigen::VectorXd build_derived_state(Eigen::VectorXd state, double t);

    void observe(Eigen::VectorXd &state , double t );

    // Function that writes states to output csv
    void write_output(std::vector<double>& time, std::vector<Eigen::VectorXd>& states, std::vector<Eigen::VectorXd>& derived_states);

    // Function that writes callback from NASA Horizons API
    static size_t write_callback(void* contents, size_t size, size_t nmemb, std::string* userp);

    // Function that builds input for nrlmsise-00 model from current state
    std::vector<double> build_atmosphere_state(Eigen::VectorXd state, double r, double t);

    // Function that calculates lla for the current state
    std::vector<double> lla(Eigen::VectorXd state, double r, double t);

    // Function that converts unix time to specific datestring for NASA Horizons API
    std::string unix_to_string(int unixTimeSeconds);

};