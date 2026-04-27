
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>


#include <json/json.h>


#include <Optimizers/GradientDescent.h>
#include <Functions.h>
#include <Simulator.h>

GradientDescent::GradientDescent(double tf, double t_burn, Json::Value spacecraft, Json::Value central_body) : 
Optimizer(tf, t_burn, spacecraft, central_body) {}

Eigen::Vector3d GradientDescent::optimize(){

    double bx_target = spacecraft["target_b_plane_coordinates"]["x"].asDouble();
    double by_target = spacecraft["target_b_plane_coordinates"]["y"].asDouble();

    Eigen::Vector3d a(0,0,0);
    return a;

}