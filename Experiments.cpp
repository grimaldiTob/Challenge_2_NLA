/*
I guess the name "Experiments" will be carried over the semester for whatever 
we're doing together.
*/

#include <iostream>
#include <tuple>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;
using MyTuple = tuple<int, int>; // define tuple type that will be used in order to store matrices...

/* Function that creates an adjacency matrix from any array of tuples that is passed to it.
    The function uses Triplets in order to store the unitary values correctly in a SparseMatrix.
    
    (Appunto: Non so se l'utilizzo della matrice sparsa sia il più adatto, visto che 
    teoricamente in un grafo tutti i nodi possono essere connessi con tutti e in quel caso non sarebbe
    più sparsa la matrice di adiacenza, però visto che in questo caso lavoriamo con cluster di dati,
    penso che comunque ci sia una sostanziale correlazione tra i dati. Ditemi che ne pensate e nel caso cambiate implementazione.)*/
SparseMatrix<double> build_adjacency_matrix(const vector<MyTuple> & arr){
    
    // we need to know the nth-biggest node that is connected to other nodes in order to initialize the sparse matrix
    int max_index = 0;
    for (const auto& tupl : arr)
        max_index = std::max({max_index, std::get<0>(tupl), std::get<1>(tupl)});
    SparseMatrix<double> mat(max_index, max_index);

    vector<Triplet<int>> tripletList;

    // of course I reserve the arr.size() * 2 because the matrix that will be returned is symmetric
    tripletList.reserve(arr.size()*2);

    for (const auto& [u, v] : arr) {
        tripletList.emplace_back(u-1, v-1, 1.0); // -1 is because I guess node 1 is in index 0 of the matrix, otherwise we have one useless row/column in the matrix
        tripletList.emplace_back(v-1, u-1, 1.0); // I dont know if this is the most efficient way but any Adjacency Matrix is symmetric
    }
    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;
}

int main(int argc, char** argv){
    using SpMat = Eigen::SparseMatrix<double>;
    using SpVec = Eigen::VectorXd;
    
    /*if(argc < 2){
        cerr << "Wrong command format!\n" << endl:
        cerr << "Correct format: ./experiments whatever it is needed I still don't know" << endl;
    } */
    
    vector<MyTuple> graph_1 = {
        MyTuple{1, 2},
        MyTuple{1, 4},
        MyTuple{2, 3},
        MyTuple{3, 4},
        MyTuple{3, 5},
        MyTuple{5, 6},
        MyTuple{6, 7},
        MyTuple{7, 9},
        MyTuple{7, 8},
        MyTuple{5, 9},
        MyTuple{8, 9},
    };

    SpMat Ag = build_adjacency_matrix(graph_1);

    // ====================================== REQUEST 1 ======================================

    cout << "Adjacency matrix size: " << Ag.rows() << "x" << Ag.cols() << endl;
    cout << "Nonzero entries: " << Ag.nonZeros() << endl;
    cout << "Frobenius norm: " << Ag.norm() << endl;

    // ====================================== REQUEST 2 ======================================

    SpVec ones = SpVec::Ones(Ag.cols());
    SpVec vg = Ag * ones; // Eigen has a method to compute the sum of the elements in a row but just in a Dense Matrix.
    // multiplying Ag for a vector of ones gives us the sum over the rows...guess we also need vector vg for later. 
    SpMat Dg = vg.asDiagonal().toDenseMatrix().sparseView();
    SpMat Lg = Dg - Ag;
    SpVec yg = Lg * vg;
    
    cout << "Norm of vector y is: " << yg.norm() << endl;

    // code copied from last challenge (check for SPD)
    Eigen::SimplicialLLT<SparseMatrix<double>> chol(Lg);
    cout << "Lg is SPD?: " << (chol.info() == Success ? "Yes" : "No") << endl; // Lg is for sure symmetric, no need to check on that 




}