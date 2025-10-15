/*
I guess the name "Experiments" will be carried over the semester for whatever 
we're doing together.
*/

#include <iostream>
#include <tuple>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Eigenvalues> 

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

using namespace std;
using namespace Eigen;
using MyTuple = tuple<int, int>; // define tuple type that will be used in order to store matrices...

// Function to print the adjacency matrix (dense format) USED FOR DEBUGGING
void print_adjacency_matrix(const Eigen::SparseMatrix<double>& mat) {
    Eigen::MatrixXd dense = Eigen::MatrixXd(mat);
    std::cout << "Adjacency Matrix:" << std::endl;
    for (int i = 0; i < dense.rows(); ++i) {
        for (int j = 0; j < dense.cols(); ++j) {
            std::cout << dense(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

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

// used for debugging 
void print_vector(const VectorXd &v) {
    cout << "[ ";
    for(int i = 0; i < v.size(); i++) {
        cout << " "<< v.coeff(i) << " ";
    }
    cout << " ]" << endl;
}

int main(int argc, char** argv){
    using SpMat = Eigen::SparseMatrix<double>;
    using SpVec = Eigen::VectorXd; // should change the name because is defined a dense vector
    
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

    print_adjacency_matrix(Ag);
    cout << "\n\n" << endl;

    // ====================================== REQUEST 1 ======================================

    
    cout << "Adjacency matrix size: " << Ag.rows() << "x" << Ag.cols() << endl;
    cout << "Nonzero entries: " << Ag.nonZeros() << endl;
    cout << "Frobenius norm: " << Ag.norm() << endl; //TODO: Ag.norm() is actually the Froebenius norm?

    // ====================================== REQUEST 2 ======================================

    SpVec ones = SpVec::Ones(Ag.cols());
    SpVec vg = Ag * ones; // Eigen has a method to compute the sum of the elements in a row but just in a Dense Matrix.
    // multiplying Ag for a vector of ones gives us the sum over the rows...guess we also need vector vg for later. 
    SpMat Dg = vg.asDiagonal().toDenseMatrix().sparseView();
    SpMat Lg = Dg - Ag;
    // SpVec yg = Lg * vg; - The assignment specifies multiplying by a vector of ones
    SpVec yg = Lg * ones; // return a zeros vector 

    cout << "yg = ";
    print_vector(yg);

    cout << "Norm of vector y is: " << yg.norm() << endl;

    // code copied from last challenge (check for SPD)
    Eigen::SimplicialLLT<SparseMatrix<double>> chol(Lg);
    cout << "Lg is SPD?: " << (chol.info() == Success ? "Yes" : "No") << endl; // Lg is for sure symmetric, no need to check on that 

    // ====================================== REQUEST 3 ======================================
    Eigen::SelfAdjointEigenSolver<SpMat> eigensolver(Lg);
    if (eigensolver.info() == Success) {
        VectorXd evals = eigensolver.eigenvalues(); // vector of eigenvalues
        MatrixXd evecs = eigensolver.eigenvectors(); // matrix with the eigenvectors as columns
        for (int i = 0; i < evals.size(); ++i) {
            if (std::abs(evals[i]) < 1.0e-15) { // clamping to zero the very small eigenvalues
                evals[i] = 0.0;
            }
        }

        cout << "\nMinimum eigenvalue is: " << evals.minCoeff() << endl;
        cout << "Maximum eigenvalue is: " << evals.maxCoeff() << "\n" << endl;

    // ====================================== REQUEST 4 ======================================
        double second_smallest = std::numeric_limits<double>::max(); // set a variable to infinity as a reference
        int index = -1;
        for (int i = 0; i < evals.size(); ++i) {
            if (evals[i] > 0 && evals[i] < second_smallest) { 
                second_smallest = evals[i];
                index = i;
            }
        }
        if(index != -1){
            cout << "The second smallest eigenvalue is: " << second_smallest << endl;
            cout << "Its corresponding eigenvector is: \n" << evecs.col(index) << endl;

            // Identifying the clusters by looking at the signs of the components of the eigenvector
            cout << "The clusters are: " << endl;
            cout << "Cluster 1 (positive components): ";
            for (int i = 0; i < evecs.rows(); ++i) // iterate over the the index column
                if (evecs(i, index) > 0)
                    cout << (i + 1) << " "; // +1 to match the node numbering
            cout << endl;
            cout << "Cluster 2 (negative components): ";
            for (int i = 0; i < evecs.rows(); ++i)
                if (evecs(i, index) < 0)
                    cout << (i + 1) << " "; // +1 to match the node numbering
            cout << endl;
        }else{
            cout << "There is no second smallest eigenvalue!" << endl;
        }
    }else{
        cerr << "Eigenvalue decomposition failed!" << endl;
    }

    // =================================== REQUEST 5 ============================================

    SpMat As;

    loadMarket(As, "social.mtx");
    cout << "\n\n\n Matrix As: \n" << endl;
    cout << "As: Adjacency matrix size =  " << As.rows() << "x" << As.cols() << endl;
    cout << "As: Nonzero entries =  " << As.nonZeros() << endl;
    cout << "As: Frobenius norm = " << As.norm() << endl;

    // ================================== REQUEST 6 =============================================

    // same routine applied in Request 2 but this time on a new adjacency matrix
    ones = SpVec::Ones(As.cols());
    SpVec vs = As * ones; 
    SpMat Ds = vs.asDiagonal().toDenseMatrix().sparseView();
    SpMat Ls = Ds - As;
    SpVec ys = Ls * ones;  

    // I need this for later
    int r = Ls.rows(); 
    int c = Ls.cols();
    int nnz = Ls.nonZeros();

    cout << "Norm of vector ys is: " << ys.norm() << endl;

    SimplicialLLT<SparseMatrix<double>> chols(Ls);
    cout << "Ls is SPD?: " << (chols.info() == Success ? "Yes" : "No") << endl;  
    cout << "Ls Nonzeros entries = " << nnz << endl; // 351 (diagonal) + 8802 (As.nonZeros()) = 9153 (Ls.nonZeros())

    // ================================== REQUEST 7 =============================================

    Ls.coeffRef(0, 0) += 0.2; // add the perturbation
    
    // saving routine, same code as the first challenge
    FILE* outLs = fopen("Ls.mtx", "w");
    fprintf(outLs,"%%%%MatrixMarket matrix coordinate real general\n");
    fprintf(outLs,"%d %d %d\n", r, c, nnz);
    for (int k=0; k<Ls.outerSize(); ++k)
    {
        for (SpMat::InnerIterator it(Ls, k); it; ++it)
        {
            fprintf(outLs, "%d %d %f\n", static_cast<int>(it.row()) + 1, static_cast<int>(it.col()) + 1, it.value()); 
        }
    }
    fclose(outLs);
    // check README.md for LIS output.
    
    return 0;
}