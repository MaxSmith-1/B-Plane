#include <iostream>
#include <json/json.h>

#include <fstream>
#include <string>

#include <Eigen/Dense>
#include <vector>




class Simulator
{

    public:

    // Initialze Constructor
    Simulator(double tf, Json::Value spacecraft, Json::Value body);
    
    // Initialize simulation loop function
    void simulate();

    private:

    // Input variables
    double tf;
    Json::Value spacecraft;
    Json::Value body;

    double abs_tol;
    double rel_tol;

    // Gravitational Constant
    double G = 6.674e-20;
    double mu;
    double r;

    // Initialize state vectors
    Eigen::VectorXd state;
    //Eigen::VectorXd derived_state;

    // Initialize time column
    std::vector<double> time;

    // Initialize list for state vector storage
    std::vector<Eigen::VectorXd> states;
    std::vector<Eigen::VectorXd> derived_states;

    // Function that gets called on every simulation loop
    void ode_function(const Eigen::VectorXd &x, Eigen::VectorXd &dxdt, const double t);

    // Function that calculates derived state values on each simulation loop
    Eigen::VectorXd build_derived_state(Eigen::VectorXd state);

    void observe(const Eigen::VectorXd &state , double t );

    // Function that writes states to output csv
    void write_output(std::vector<double>& time, std::vector<Eigen::VectorXd>& states, std::vector<Eigen::VectorXd>& derived_states);

};