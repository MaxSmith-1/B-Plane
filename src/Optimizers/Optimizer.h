#include <iostream>
#include <Eigen/Dense>
#include <json/json.h>
#include <./Simulator.h>

// Class that inputs a simulator class with a satellite config, outputs a burn to insert close to target impact parameter at t_burn

class Optimizer{

    public:

    Optimizer(double tf, double t_burn, Json::Value spacecraft, Json::Value central_body) 
        : t_burn(t_burn), tf(tf), spacecraft(spacecraft), central_body(central_body),
        sim(tf, this->spacecraft, central_body, false)  // ← this->spacecraft is the member
    {
        Eigen::Vector2d tbc(
            this->spacecraft["target_b_plane_coordinates"]["X"].asDouble(), 
            this->spacecraft["target_b_plane_coordinates"]["Y"].asDouble()
        );
        target_b_coordinates = tbc;
    }

    virtual Eigen::Vector3d optimize() = 0;

    virtual ~Optimizer() {}

    protected:
        double tf;
        double t_burn;
        Json::Value spacecraft;   // ← must come before sim
        Json::Value central_body;
        Simulator sim;            // ← declared after spacecraft
        Eigen::Vector2d target_b_coordinates;
        unsigned int i;





};
