#include <iostream>
#include <json/json.h>
#include <Optimizers/Optimizer.h>
#include <Eigen/Dense>


class GradientDescent : public Optimizer {

    public:

    GradientDescent(double tf, double t_burn, Json::Value spacecraft, Json::Value central_body);

    Eigen::Vector3d optimize() override;


};