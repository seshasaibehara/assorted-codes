#include "headers.h"
#include "../eigen-git-mirror/Eigen/Core"
#include <iostream>
#include <fstream>
#include <chrono>

int main(int argc, char* argv[]){
    input input_data;
    input_data.point = 0;
    input_data.pair = 0.04;
    input_data.size_x = 200;
    input_data.size_y = 200;
    
    auto start = std::chrono::high_resolution_clock::now();
    std::ofstream outfile;
    outfile.open("output.txt");
    std::cout << "Job started" << std::endl;
    Eigen::MatrixXi initial_mircrostate = generate_initial_microstate(input_data.size_x, input_data.size_y);
    for (double mu=-0.5; mu<=0.5; mu=mu+0.025){
        for (int T=300; T<2000; T=T+50){
            std::tuple<Eigen::MatrixXi,double,double,double> final_microstate = metropolis_algorithm(input_data,initial_mircrostate,mp,mu,T);
            initial_mircrostate = std::get<0>(final_microstate); 
        //writing data out to a file
            outfile << mu << "  " << T << " " << std::get<1>(final_microstate) << " "<< std::get<2>(final_microstate) << "  "<< std::get<3>(final_microstate) << std::endl;
          std::cout << "Finished run for temperature " << T << " and chemical potential " << mu << std::endl;         
  
        }
    }
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    std::cout<< duration.count() << std::endl;

    outfile.close();

    return 0;
}
