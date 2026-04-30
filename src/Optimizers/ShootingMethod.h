#include <iostream>
#include <json/json.h>
#include <Optimizers/Optimizer.h>
#include <Eigen/Dense>


class ShootingMethod : public Optimizer {

    public:

    ShootingMethod(double tf, double t_burn, Json::Value spacecraft, Json::Value central_body);

    Eigen::Vector3d optimize() override;


    private:

    void add_burn(Eigen::Vector3d velocity);

    std::array<double, 2> get_b_coordinates(double tf, Json::Value spacecraft, Json::Value central_body, bool monte_carlo);

    Eigen::Vector3d get_rtn_burn(Eigen::Vector3d velocity, double tf, Json::Value spacecraft, Json::Value central_body, bool monte_carlo);

    Eigen::Matrix<double, 2, 3> get_Jacobean(Eigen::Vector3d velocity);

    Eigen::Vector2d get_residual(Eigen::Vector3d velocitym, double bx_target, double by_target);




};