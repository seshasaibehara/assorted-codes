#include "../eigen-git-mirror/Eigen/Dense"
#include "headers.h"
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <cmath>

Vector3d cart2fract(Vector3d cart_coords, Matrix3d lattice){
    Vector3d fract_coords;
    fract_coords = lattice.inverse()*cart_coords;
    return fract_coords;
}
    
Vector3d fract2cart(Vector3d fract_coords, Matrix3d lattice){
    Vector3d cart_coords;
    cart_coords = lattice*fract_coords;
    return cart_coords;
}

bool _is_coordinate_inside_crystal(double coordinate){
    if (coordinate<-PREC || coordinate>=1-PREC){
        return false;
    }
    return true;
}

void move_the_atom_inside_the_crystal(Vector3d* fract_coords_ptr, Matrix3d lattice){
    auto& fract_coords = *fract_coords_ptr;
    for (int i=0; i<3; ++i){
        if (_is_coordinate_inside_crystal(fract_coords(i))==false){
            fract_coords(i) = fract_coords(i)-floor(fract_coords(i)+PREC);        
        }
    }
    return;
}

bool compare_doubles(double a, double b, double absEpsilon, double relEpsilon){
    // Check if the numbers are really close -- needed when comparing numbers near zero.
    double diff{fabs(a - b)};
    if (diff <= absEpsilon)
        return true;
    // Otherwise fall back to Knuth's algorithm
    return diff <= ((fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * relEpsilon);
    //
}


crystal read_input_data(std::string file_name){
    
    std::ifstream input;
    input.open(file_name, std::ios::in);
    std::string line;

    crystal crystal_data;
    Matrix3d lattice_data;
    basis_class temp_basis;
    Vector3d temp_vector;
    std::vector<basis_class> basis_data; 
    
    //Line 1
    getline(input, line);
    std::istringstream ss0(line);
    double x0,y0,z0;
    ss0 >> x0 >> y0 >> z0;
    lattice_data(0,0) = x0;
    lattice_data(1,0) = y0;
    lattice_data(2,0) = z0;
    //Line 2
    getline(input, line);
    std::istringstream ss1(line);
    double x1,y1,z1;
    ss1 >> x1 >> y1 >> z1;
    lattice_data(0,1) = x1;
    lattice_data(1,1) = y1;
    lattice_data(2,1) = z1;
    //Line 3
    getline(input, line);
    std::istringstream ss2(line);
    double x2,y2,z2;
    ss2 >> x2 >> y2 >> z2;
    lattice_data(0,2) = x2;
    lattice_data(1,2) = y2;
    lattice_data(2,2) = z2;
    //From line 4 onwards, the format is basis so use while 
    while (getline(input, line)){
        std::istringstream ss(line);
        std::string s;
        double x,y,z;
        
        ss >> s >> x >> y >> z;
        
        temp_basis.atom_type = s;
        temp_vector(0) = x;
        temp_vector(1) = y;
        temp_vector(2) = z;
        temp_basis.coordinates = temp_vector;

        basis_data.push_back(temp_basis);
    }

    input.close();

    crystal_data.lattice = lattice_data;
    crystal_data.basis = basis_data;

    //Printing to check whether it read correctly
//    std::cout << crystal_data.lattice << std::endl;
//    for (auto &x: crystal_data.basis){
//        std::cout << x.atom_type << std::endl;
//        std::cout << x.coordinates << std::endl;
//    }

    return crystal_data;
}

std::vector<Vector3d> calc_grid_points(const int radius, const Matrix3d lattice_data){

    std::vector<Vector3d> grid_points;
    Vector3d grid_point;
    for (int i=-radius; i<=radius; ++i){
        for (int j=-radius; j<=radius; ++j){
            for (int k=-radius; k<=radius; ++k){
                grid_point << i*lattice_data(0,0) + j*lattice_data(0,1) + k*lattice_data(0,2),
                    i*lattice_data(1,0) + j*lattice_data(1,1) + k*lattice_data(1,2),
                    i*lattice_data(2,0) + j*lattice_data(2,1) + k*lattice_data(2,2);
                grid_points.push_back(grid_point);
            }
        }
    }
    return grid_points;
}

std::vector<Matrix3d> calc_l_primes(const std::vector<Vector3d>& grid_points){
    std::vector<Matrix3d> lprimes;
    Matrix3d lprime;

    for (auto &x: grid_points){
        for (auto &y: grid_points){
            for (auto &z: grid_points){
                lprime << x(0), y(0), z(0),
                        x(1), y(1), z(1),
                        x(2), y(2), z(2);
                lprimes.push_back(lprime);
            }
        }
    }
    return lprimes;
}

bool is_sym_op(const Matrix3d possible_sym_op){
    if(compare_doubles(possible_sym_op.determinant(),1,PREC,PREC) == true || compare_doubles(possible_sym_op.determinant(),-1,PREC,PREC)==true){
        Matrix3d check_identity;
        check_identity = possible_sym_op.transpose()*possible_sym_op;
        if(check_identity.isIdentity(PREC)==true){
            return true;
        }
    }
    return false;
}

int count_positive_unit_eigen_values(sym_op_class sym_op_in) {

    eigensolver ev(sym_op_in.coord_matrix, false);
    Vector3d evs;
    evs = ev.eigenvalues().real();
    int counter = 0;
    for (int i=0; i<3; i++) {
        if (compare_doubles(evs(i,0),1,PREC,PREC)){
        counter = counter + 1;
        }      
    }
    return counter;
}

//Generate labels for various sym ops based on the conditions
std::string label_sym_ops(sym_op_class sym_op_in){
        
    double det = sym_op_in.coord_matrix.determinant();
    double trace = sym_op_in.coord_matrix.trace();
    int eigen_values = count_positive_unit_eigen_values(sym_op_in);
    std::string label;

    if (compare_doubles(det,1,PREC,PREC)==true && compare_doubles(trace,3,PREC,PREC)==true){
        label = "Identity";    
    }
        
    if (compare_doubles(det,-1,PREC,PREC)==true && compare_doubles(trace,-3,PREC,PREC)==true){
        label = "Inverse";
    }

    if (compare_doubles(det,1,PREC,PREC)==true && eigen_values==1){
        label = "Rotation";
    }

    if (compare_doubles(det,-1,PREC,PREC)==true && compare_doubles(trace,-3,PREC,PREC)==false && eigen_values==0){
        label = "Improper Rotation";
    }

    if (compare_doubles(det,-1,PREC,PREC)==true && eigen_values==2){
        label = "Mirror";
    }
    return label;
}
//write code for printing the multiplication table of the group 
MatrixXi generate_mult_table(const std::vector<Matrix3d>& point_group){
    const int order = point_group.size();
    MatrixXi mult_table(order, order);
    
    for (int i = 0; i < order; i++){
        for (int j = 0; j < order; j++){
            const Matrix3d product = point_group[i]*point_group[j];
            auto it_k = std::find_if(point_group.begin(), point_group.end(),
                                     [&](const Matrix3d& symop) { return symop.isApprox(product, PREC); });
            if (it_k != point_group.end()){
                mult_table(i,j) = it_k - point_group.begin();
            }
            else{
                std::cerr << "Group is not closed" << std::endl;
                return MatrixXi::Constant(order, order, -1);
            }
        }
    }
    return mult_table;
}

std::vector<sym_op_class> calc_point_group_of_lattice(const crystal crystal_data){
    
    std::vector<sym_op_class> point_group;
    sym_op_class temp_sym_op;
    std::vector<Vector3d> grid_points = calc_grid_points(1,crystal_data.lattice);
    std::vector<Matrix3d> lprimes = calc_l_primes(grid_points);
    Matrix3d possible_sym_op;
    std::vector<Matrix3d> syms;    

    for (auto &x: lprimes){
        possible_sym_op = x*crystal_data.lattice.inverse();
        if (is_sym_op(possible_sym_op)==true){
            temp_sym_op.coord_matrix = possible_sym_op;
            temp_sym_op.label = label_sym_ops(temp_sym_op);
            point_group.push_back(temp_sym_op);
            syms.push_back(possible_sym_op);
            //Printing sym ops
            //std::cout << temp_sym_op.label << std::endl;
            //std::cout << possible_sym_op << std::endl;
            //std::cout << "----" << std::endl;
        }
    }
    //std::cout << symmetry_operations.size() << std::endl;
    MatrixXi mult_table = generate_mult_table(syms);
    std::cout << "Multiplication table of point group:" << std::endl;
    std::cout << mult_table << std::endl;
    std::cout << "---------------" << std::endl;
    return point_group;
        
}


bool are_vectors_equal(Vector3d v0, Vector3d v1){
    if (compare_doubles(v0(0),v1(0),PREC,PREC)==true && compare_doubles(v0(1),v1(1),PREC,PREC)==true && compare_doubles(v0(2),v1(2),PREC,PREC)){
        return true;
    }

    return false;
}

bool find_if(const std::vector<Vector3d>& possible_translations, Vector3d possible_translation){
    for (auto &x:possible_translations){
        if (are_vectors_equal(x,possible_translation)){
                return true;
        }
    }
    return false;
}

bool is_structure_equivalent(const std::vector<basis_class>& shifted_basis, const std::vector<basis_class>& old_basis){
    
    for (auto &x: shifted_basis){
        int counter = 0;
        for (auto &y: old_basis){
            if (x.atom_type==y.atom_type && are_vectors_equal(x.coordinates,y.coordinates)==true){
                counter = counter + 1;
            }
        }
        if (counter==0){
            return false;
        }
    }
    return true;
}

std::vector<basis_class> generate_shifted_basis(const std::vector<basis_class>& new_basis, Vector3d possible_translation, Matrix3d lattice){
    std::vector<basis_class> shifted_basis;
    basis_class temp_basis;
//    std::cout << possible_translation << std::endl;
    for (auto &x: new_basis){
        temp_basis.atom_type = x.atom_type;
        temp_basis.coordinates = x.coordinates + possible_translation;//        
//        std::cout << temp_basis.coordinates << std::endl;
// Add code if the point falls outside the box, move it to inside the box
        temp_basis.coordinates = cart2fract(temp_basis.coordinates, lattice);
        move_the_atom_inside_the_crystal(&temp_basis.coordinates, lattice); 
        temp_basis.coordinates = fract2cart(temp_basis.coordinates, lattice);
        //        std::cout << temp_basis.coordinates << std::endl;
        shifted_basis.push_back(temp_basis);
    }
    return shifted_basis;
}

std::vector<Vector3d> generate_translations(const std::vector<basis_class>& new_basis, const std::vector<basis_class>& old_basis, const Matrix3d lattice){
    std::vector<Vector3d> translations;
    std::vector<Vector3d> possible_translations;
    std::vector<basis_class> shifted_basis;
    Vector3d temp_vector;
    
        for (auto &x: new_basis){
            temp_vector = old_basis[0].coordinates - x.coordinates;
            temp_vector = cart2fract(temp_vector, lattice);
            move_the_atom_inside_the_crystal(&temp_vector,lattice);            
            //std::cout << temp_vector << std::endl;
            if (find_if(possible_translations,temp_vector)==false){
                possible_translations.push_back(temp_vector);
            }
    }
    for (auto &possible_translation: possible_translations){
        
        std::cout << possible_translation << std::endl;
        shifted_basis = generate_shifted_basis(new_basis, fract2cart(possible_translation,lattice), lattice);
        std::cout << is_structure_equivalent(shifted_basis,old_basis) << std::endl;
        std::cout << "----" << std::endl;
        if (is_structure_equivalent(shifted_basis, old_basis)==true){
            translations.push_back(possible_translation);
        }
    }
    return translations;
}

std::vector<sym_op_class> calc_factor_group(const std::vector<sym_op_class>& point_group, const crystal crystal_data){
    

    std::vector<sym_op_class> factor_group;
    sym_op_class temp_sym_op;
    std::vector<basis_class> new_basis;
    basis_class temp_basis;

    std::vector<Vector3d> translations;

    for (auto &symmetry_operation: point_group){
        for (auto &basis: crystal_data.basis){
            temp_basis.coordinates = symmetry_operation.coord_matrix*basis.coordinates;
            temp_basis.atom_type = basis.atom_type;
            //std::cout << temp_basis.coordinates << std::endl;
            new_basis.push_back(temp_basis);  
         }

        std::cout << symmetry_operation.coord_matrix << std::endl; 
        std::cout << symmetry_operation.label << std::endl;
        translations = generate_translations(new_basis,crystal_data.basis,crystal_data.lattice);
        std::cout << "====" << std::endl;
        for (auto &translation: translations){
            temp_sym_op.coord_matrix = symmetry_operation.coord_matrix;
            temp_sym_op.translation = translation;
            temp_sym_op.label = symmetry_operation.label;
            factor_group.push_back(temp_sym_op);

        }
        new_basis.clear();
        new_basis.resize(0);
    }
    return factor_group;
}

void print_point_group(const std::vector<sym_op_class>& point_group){
    
    for (auto &x: point_group){
        std::cout << x.coord_matrix.format(CleanFmt) << std::endl;
        std::cout << x.label << std::endl;
        
        std::cout << "-----------" << std::endl;
    }
     
    std::cout << "Point group operations: " << point_group.size() << std::endl;

    std::cout << "============================" << std::endl;
}

void print_factor_group(const std::vector<sym_op_class>& factor_group){
   
    for (auto &x: factor_group){
        std::cout << x.label << std::endl; 
        std::cout << "S: " << std::endl;
        std::cout << x.coord_matrix.format(CleanFmt) << std::endl;
        std::cout << "\u03A4: " << std::endl;
        std::cout << x.translation.format(CleanFmt) << std::endl;
        
        std::cout << "---------" << std::endl;
    }
     
    std::cout << "Factor group operations: " << factor_group.size() << std::endl;
    std::cout << "==================================" << std::endl;
}


int main(int argc, char* argv[]){

    if (argc == 1){
        std::cerr << "No crystal specified" << std::endl;
        return 1;
    }

    crystal crystal_data = read_input_data(argv[1]);
    std::vector<sym_op_class> point_group = calc_point_group_of_lattice(crystal_data);
    std::vector<sym_op_class> factor_group = calc_factor_group(point_group, crystal_data);
//    print_point_group(point_group);
    print_factor_group(factor_group);

    return 0;

}
