#include "../eigen-git-mirror/Eigen/Dense"
#include "headers.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <numeric>
#define kb 8.6e-5;

//input read_input_file(std::string file_name){
//    std::ifstream in;
//    in.open(file_name, std::ios::in);
//    std::string line;
//    input input_data;
//
//    getline(in,line);
//    std::istringstream ss(line);
//    double x,y;
//    ss >> x >> y;
//    input_data.point = x;
//    input_data.pair = y;
//    
//    getline(in,line);
//    std::istringstream ss1(line);
//    int Nx, Ny;
//    ss1 >> Nx >> Ny;
//    input_data.size_x = Nx;
//    input_data.size_y = Ny;
//    
//    std::vector<int> test;
//    int p;
//    getline(in,line);
//    std::istringstream sst(line);
//    while (sst >> p){
//        test.push_back(p);
//    }
//   
//    for (auto &x: test){
//        std::cout << x << std::endl;
//    }
//    return input_data;
//}


Eigen::MatrixXi generate_initial_microstate(int size_x, int size_y){
   Eigen::MatrixXi initial_mircrostate = Eigen::MatrixXi::Constant(size_x, size_y, -1); 
   return initial_mircrostate;
}

std::tuple<Eigen::MatrixXi,int,int> generate_perturbation(const Eigen::MatrixXi& previous_microstate){
    Eigen::MatrixXi perturbed_microstate = previous_microstate;
    int i = std::rand()/((RAND_MAX + 1u)/previous_microstate.rows());
    int j = std::rand()/((RAND_MAX + 1u)/previous_microstate.cols());
    perturbed_microstate(i,j) = -previous_microstate(i,j);
    std::tuple<Eigen::MatrixXi,int,int> current_state = std::make_tuple(perturbed_microstate,i,j);
    return current_state;
}

int mod(const Eigen::MatrixXi& state, int i, int j, int Nx, int Ny){
    if (i==Nx){
        i = 0;
    }
    if (j==Ny){
        j = 0;
    }
    if (i==-1){
        i = Nx-1;
    }
    if (j==-1){
        j = Ny-1;
    }
    return state(i,j);
}

double compute_deltaE(const Eigen::MatrixXi& previous_microstate,const std::tuple<Eigen::MatrixXi,int,int>& current_state, double Vp, double Vnn){
    int i = std::get<1>(current_state);
    int j = std::get<2>(current_state);
    int Nx = previous_microstate.rows();
    int Ny = previous_microstate.cols();
    double deltaE = -2*previous_microstate(i,j)*(Vp+Vnn*(mod(previous_microstate,i+1,j,Nx,Ny))+Vnn*(mod(previous_microstate,i,j+1,Nx,Ny))+Vnn*(mod(previous_microstate,i-1,j,Nx,Ny))+Vnn*(mod(previous_microstate,i,j-1,Nx,Ny)));
//    std::cout << deltaE << std::endl;
    return deltaE;
}

double compute_deltaN(const Eigen::MatrixXi& previous_microstate,const std::tuple<Eigen::MatrixXi,int,int>& current_state){
    int i = std::get<1>(current_state);
    int j = std::get<2>(current_state);
    double deltaN = -previous_microstate(i,j);
    return deltaN;
}

double compute_deltaOmega(double deltaE, double deltaN, double mu){
    double deltaOmega = deltaE - (mu*deltaN);
//    std::cout << deltaOmega << std::endl;
    return deltaOmega;
}

bool is_acceptable_microstate(const Eigen::MatrixXi& previous_microstate,const Eigen::MatrixXi& current_state, double deltaOmega, double T){
    if (deltaOmega < 0){
        return true;
    }
    double beta = 1/kb;
    beta = beta/T;
    double prob = std::exp(-beta*deltaOmega);
    double random_number = std::rand()/((double)RAND_MAX + 1u);
//    std::cout << prob << " " << random_number << std::endl;
    if ( prob > random_number){
        return true;
    }
    return false;
}
double compute_energy(const Eigen::MatrixXi& microstate,const input input_data){
    double energy = 0;
    for (int i=0; i<microstate.rows(); ++i){
        for (int j=0; j<microstate.cols(); ++j){
            energy = energy + input_data.pair*microstate(i,j)*(mod(microstate,i-1,j,input_data.size_x,input_data.size_y) + mod(microstate,i,j+1,input_data.size_x,input_data.size_y));
        }
    }
//    std::cout << energy << std::endl;
    return energy;
}
double compute_composition(const Eigen::MatrixXi& microstate){
    double composition = 0;
    for (int i=0; i<microstate.rows(); ++i){
        for (int j=0; j<microstate.cols(); ++j){
            if (microstate(i,j)==1){
                composition = composition + 1;
            }
        }
    }
    return composition;
}

double add_elements_in_vector(const std::vector<double>& v){
    double sum = 0;
    for (auto &x: v){
        sum = sum + x;
    }
    return sum;
}

std::tuple<Eigen::MatrixXi,double,double,double> metropolis_algorithm(const input input_data,Eigen::MatrixXi& microstate,int mp, double mu, double T) {
    std::tuple<Eigen::MatrixXi,int,int> perturbed_microstate; 
    std::vector<double> energies;
    std::vector<double> compositions;
    std::vector<double> q1;
    std::vector<double> q2;
    int size = input_data.size_x*input_data.size_y;
    int i,j=0;
    while (j<mp){
        while(i<size){
            perturbed_microstate = generate_perturbation(microstate);
//            std::cout << std::get<0>(perturbed_microstate) << std::endl;
            double deltaE = compute_deltaE(microstate,perturbed_microstate,input_data.point,input_data.pair);
//            std::cout << deltaE << std::endl;
            double deltaN = compute_deltaN(microstate,perturbed_microstate);
            double deltaOmega = compute_deltaOmega(deltaE,deltaN,mu);
//            std::cout << deltaOmega << std::endl;
//            std::cout << is_acceptable_microstate(microstate,std::get<0>(perturbed_microstate),deltaOmega,T) << std::endl;
//            std::cout << "-----" << std::endl;
            if (is_acceptable_microstate(microstate,std::get<0>(perturbed_microstate),deltaOmega,T)==true){
                microstate = std::get<0>(perturbed_microstate);
            }
                        ++i;
        }
        if(j>2000){    
                double energy = compute_energy(microstate,input_data);
                double composition = compute_composition(microstate);
                energies.push_back(energy);
                compositions.push_back(composition/(double)(input_data.size_x*input_data.size_y));
                q1.push_back((energy-(mu*composition))*(energy-(mu*composition)));
                q2.push_back((energy-(mu*composition)));
            }
    ++j;
    }
    double average_energy = add_elements_in_vector(energies)/(double)energies.size();
    double average_composition = add_elements_in_vector(compositions)/(double)compositions.size();
    double heat_capacity = add_elements_in_vector(q1)/q1.size();
    double q3 = add_elements_in_vector(q2)/q2.size();
    heat_capacity = (heat_capacity - (q3*q3))/((8.6e-5)*T*T);
    std::tuple<Eigen::MatrixXi,double,double,double> results_to_return = std::make_tuple(microstate,average_energy,average_composition,heat_capacity);
    return results_to_return; 
}


