#ifndef HEADERS_H
#define HEADERS_H
#include "../eigen-git-mirror/Eigen/Dense"

struct input{
    double point, pair;
    int size_x, size_y;
};

//    double mu, T;
//    int mp;
//    input read_input_file(std::string file_name);
    Eigen::MatrixXi generate_initial_microstate(int size_x, int size_y);
    std::tuple<Eigen::MatrixXi,int,int> generate_perturbation(Eigen::MatrixXi& previous_microstate);
    double compute_deltaE(Eigen::MatrixXi& previous_microstate, std::tuple<Eigen::MatrixXi,int,int>& current_state, double Vp, double Vnn);
    double compute_deltaN(Eigen::MatrixXi& previous_microstate, std::tuple<Eigen::MatrixXi,int,int>& current_state);
    double compute_deltaOmega(double deltaE, double deltaN, double mu);
    bool is_acceptable_microstate(Eigen::MatrixXi& previous_microstate, Eigen::MatrixXi& current_state, double deltaOmega, double T);
    std::tuple<Eigen::MatrixXi,double,double,double> metropolis_algorithm(input input_data, Eigen::MatrixXi& initial_mircrostate,int mp, double mu, double T);

    int mod(Eigen::MatrixXi& microstate, int i, int j, int Nx, int Ny); 
#endif
