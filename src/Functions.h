#include <iostream>
#include <json/json.h>

#include <fstream>
#include <string>

#include <Eigen/Dense>
#include <vector>

class Functions{

    public:

    static Eigen::Vector3d rtn_to_icrf_delta_v(Eigen::VectorXd state, const Eigen::Vector3d& delta_v_rtn){

        // Calculate r hat vector
        Eigen::Vector3d r_vec(state[0], state[1], state[2]);
        Eigen::Vector3d r_hat = r_vec / r_vec.norm();
        
        // Calculate n hat vector
        Eigen::Vector3d v_vec(state[3], state[4], state[5]);
        Eigen::Vector3d h_vec = r_vec.cross(v_vec);
        Eigen::Vector3d n_hat = h_vec / h_vec.norm();

        // Calculate t hat vector
        Eigen::Vector3d t_hat = n_hat.cross(r_vec / r_vec.norm());

        // Assemble rtn to eci rotation matrix
        Eigen::Matrix3d RTN_to_ICRF;
        RTN_to_ICRF << r_hat[0], t_hat[0],   n_hat[0],
            r_hat[1],    t_hat[1], n_hat[1],
            r_hat[2],      t_hat[2],     n_hat[2];
        

        // Return rotation product
        return RTN_to_ICRF * delta_v_rtn;

    }



    static Eigen::Vector3d icrf_to_rtn_delta_v(Eigen::VectorXd state, const Eigen::Vector3d& delta_v_icrf){

        // Calculate r hat vector
        Eigen::Vector3d r_vec(state[0], state[1], state[2]);
        Eigen::Vector3d r_hat = r_vec / r_vec.norm();
        
        // Calculate n hat vector
        Eigen::Vector3d v_vec(state[3], state[4], state[5]);
        Eigen::Vector3d h_vec = r_vec.cross(v_vec);
        Eigen::Vector3d n_hat = h_vec / h_vec.norm();

        // Calculate t hat vector
        Eigen::Vector3d t_hat = n_hat.cross(r_vec / r_vec.norm());

        // Assemble rtn to eci rotation matrix
        Eigen::Matrix3d ICRF_to_RTN;
        ICRF_to_RTN << r_hat[0], t_hat[0],   n_hat[0],
            r_hat[1],    t_hat[1], n_hat[1],
            r_hat[2],      t_hat[2],     n_hat[2];
        
        ICRF_to_RTN.transposeInPlace();

        // Return rotation product
        return ICRF_to_RTN * delta_v_icrf;

    }


};

