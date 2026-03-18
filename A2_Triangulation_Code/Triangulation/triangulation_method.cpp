/**
 * Copyright (C) 2015 by Liangliang Nan (liangliang.nan@gmail.com)
 * https://3d.bk.tudelft.nl/liangliang/
 *
 * This file is part of Easy3D. If it is useful in your research/work,
 * I would be grateful if you show your appreciation by citing it:
 * ------------------------------------------------------------------
 *      Liangliang Nan.
 *      Easy3D: a lightweight, easy-to-use, and efficient C++
 *      library for processing and rendering 3D data. 2018.
 * ------------------------------------------------------------------
 * Easy3D is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License Version 3
 * as published by the Free Software Foundation.
 *
 * Easy3D is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include "triangulation.h"
#include "matrix_algo.h"
#include <easy3d/optimizer/optimizer_lm.h>


using namespace easy3d;


/**
 * TODO: Finish this function for reconstructing 3D geometry from corresponding image points.
 * @return True on success, otherwise false. On success, the reconstructed 3D points must be written to 'points_3d'
 *      and the recovered relative pose must be written to R and t.
 */
bool Triangulation::triangulation(
        double fx, double fy,     /// input: the focal lengths (same for both cameras)
        double cx, double cy,     /// input: the principal point (same for both cameras)
        double s,                 /// input: the skew factor (same for both cameras)
        const std::vector<Vector2D> &points_0,  /// input: 2D image points in the 1st image.
        const std::vector<Vector2D> &points_1,  /// input: 2D image points in the 2nd image.
        std::vector<Vector3D> &points_3d,       /// output: reconstructed 3D points
        Matrix33 &R,   /// output: 3 by 3 matrix, which is the recovered rotation of the 2nd camera
        Vector3D &t    /// output: 3D vector, which is the recovered translation of the 2nd camera
) const
{
    /// NOTE: there might be multiple workflows for reconstructing 3D geometry from corresponding image points.
    ///       This assignment uses the commonly used one explained in our lecture.
    ///       It is advised to define a function for the sub-tasks. This way you have a clean and well-structured
    ///       implementation, which also makes testing and debugging easier. You can put your other functions above
    ///       'triangulation()'.

    std::cout << "\nTODO: implement the 'triangulation()' function in the file 'Triangulation/triangulation_method.cpp'\n\n";

    std::cout << "[Liangliang]:\n"
                 "\tSimilar to the first assignment, basic linear algebra data structures and functions are provided in\n"
                 "\tthe following files:\n"
                 "\t    - Triangulation/matrix.h: handles matrices of arbitrary dimensions and related functions.\n"
                 "\t    - Triangulation/vector.h: manages vectors of arbitrary sizes and related functions.\n"
                 "\t    - Triangulation/matrix_algo.h: contains functions for determinant, inverse, SVD, linear least-squares...\n"
                 "\tFor more details about these data structures and a complete list of related functions, please\n"
                 "\trefer to the header files mentioned above.\n\n"
                 "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
                 "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
                 "\tFor your final submission, adhere to the following guidelines:\n"
                 "\t    - submit ONLY the 'Triangulation/triangulation_method.cpp' file.\n"
                 "\t    - remove ALL unrelated test code, debugging code, and comments.\n"
                 "\t    - ensure that your code compiles and can reproduce your results WITHOUT ANY modification.\n\n" << std::flush;

    /// Below are a few examples showing some useful data structures and APIs.

    /// define a 2D vector/point
   // Vector2D b(1.1, 2.2);

    /// define a 3D vector/point
   // Vector3D a(1.1, 2.2, 3.3);

    /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
   // Vector2D p = a.cartesian();

    /// get the Homogeneous coordinates of p
   // Vector3D q = p.homogeneous();

    /// define a 3 by 3 matrix (and all elements initialized to 0.0)
   // Matrix33 A;

    /// define and initialize a 3 by 3 matrix
   // Matrix33 T(1.1, 2.2, 3.3,
               //0, 2.2, 3.3,
               //0, 0, 1);

    /// define and initialize a 3 by 4 matrix
    //Matrix34 M(1.1, 2.2, 3.3, 0,
              // 0, 2.2, 3.3, 1,
              // 0, 0, 1, 1);

    /// set first row by a vector
   // M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));

    /// set second column by a vector
   // M.set_column(1, Vector3D(5.5, 5.5, 5.5));

    /// define a 15 by 9 matrix (and all elements initialized to 0.0)
   // Matrix W(15, 9, 0.0);
    /// set the first row by a 9-dimensional vector
//W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>

    /// get the number of rows.
    //int num_rows = W.rows();

    /// get the number of columns.
   // int num_cols = W.cols();

    /// get the the element at row 1 and column 2
    //double value = W(1, 2);

    /// get the last column of a matrix
   // Vector last_column = W.get_column(W.cols() - 1);

    /// define a 3 by 3 identity matrix
    //Matrix33 I = Matrix::identity(3, 3, 1.0);

    /// matrix-vector product
   // Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4

    ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).

    if (points_0.size() < 8 || points_1.size() < 8 || points_0.size() != points_1.size()) {
        std::cout << "input is invalid" << std::endl;
        return false;
    }

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;

    //Making a copy of the points because it's a const
    std::vector<Vector2D> point_0 = points_0;
    std::vector<Vector2D> point_1 = points_1;

    // Before building matrix W we need to normalize our points:

    ////////////////////////////////////////// Normalizing Camera 1 ///////////////////////////////////////////////////
    double sum_u = 0.0, sum_v = 0.0;
    int n = point_0.size();

    // Getting the average centroid for u and v:
    for (int i = 0; i < n; i++) {
        sum_u = sum_u + point_0[i].x();
        sum_v = sum_v + point_0[i].y();
    }

    double mean_u = sum_u/n;
    double mean_v = sum_v/n;

    //To translate - we subtract the centroid from every point

    for (int i = 0; i < n; i++) {
        point_0[i].x() = point_0[i].x() - mean_u;
        point_0[i].y() = point_0[i].y() - mean_v;
    }

    //Calculating the average distance of points and scaling it to sqrt2 around the origin

    double sum_distances = 0;

    for (int i = 0; i < n; i++) {
        double u = point_0[i].x();
        double v = point_0[i].y();

        double distance = sqrt(u*u + v*v); //pythagorean distance from origin to point i
        sum_distances = sum_distances + distance;
    }

    double average_distance = sum_distances/n;
    //eg: if the current pixel is at 200, scale = sqrt2/200 = 0.00707 so we multiply every point by 0.00707 = sqrt2
    double scale = sqrt(2.0)/average_distance;

    //Scaling every point:

    for (int i = 0; i < n; i++) {
        point_0[i].x() = point_0[i].x() * scale;
        point_0[i].y() = point_0[i].y() * scale;
    }

    /////////////////////////////////////////// Normalizing Camera 2 /////////////////////////////////////////////////

    double sum_u1 = 0.0, sum_v1 = 0.0;
    int n1 = point_1.size();

    // Getting the average centroid for u and v:

    for (int i = 0; i < n1; i++) {
        sum_u1 = sum_u1 + point_1[i].x();
        sum_v1 = sum_v1 + point_1[i].y();
    }

    double mean_u1 = sum_u1/n1;
    double mean_v1 = sum_v1/n1;

    //To translate - we subtract the centroid from every point

    for (int i = 0; i<n1; i++) {
        point_1[i].x() = point_1[i].x() - mean_u1;
        point_1[i].y() = point_1[i].y() - mean_v1;
    }

    //Calculating the average distance of points and scaling it to sqrt2 around the origin

    double sum_distance1 = 0.0;

    for (int i = 0; i < n1; i++) {
        double u = point_1[i].x();
        double v = point_1[i].y();

        double distance = sqrt(u*u + v*v);
        sum_distance1 = sum_distance1 + distance;
    }

    double average_distance1 = sum_distance1/n1;
    double scale1 = sqrt(2.0)/average_distance1;

    //Scaling every point:

    for (int i = 0; i < n1; i++) {
        point_1[i].x() = point_1[i].x() * scale1;
        point_1[i].y() = point_1[i].y() * scale1;
    }

    //////////////////////////////////////////////////// Matrix W and F_q /////////////////////////////////////////////

    //Initializing matrix W:

    Matrix W(n,9,0.0);

    for (int i = 0; i < n; i++) {
        double u = point_0[i].x();
        double v = point_0[i].y();
        double u1 = point_1[i].x();
        double v1 = point_1[i].y();

        W.set_row(i,{u*u1,v*u1,u1,u*v1,v*v1,v1,u,v,1});
    }

    //Setting the SVD:
    Matrix U(n,n,0.0);
    Matrix S(n,9,0.0);
    Matrix V(9,9,0.0);

    svd_decompose(W,U,S,V);

    Vector f = V.get_column(V.cols()-1);

    Matrix33 F_q(f[0],f[1],f[2],
               f[3],f[4],f[5],
               f[6],f[7],f[8]
               );

    //////////////////////////////////////////// Constraint enforcement ///////////////////////////////////////////////

    ///Slide 48 F_q = UDVtranspose, so we do SVD on F_q so find the singular values of F_q so we can fix its rank

    Matrix U_q(3,3,0.0);
    Matrix D_q(3,3,0.0);
    Matrix V_q(3,3,0.0);

    svd_decompose(F_q,U_q,D_q,V_q);

    ///Setting d3 (2,2) to 0
    D_q(2,2) = 0.0;

    ///Recomposing the F_q : F = U_q * D_q * V_q transpose
    Matrix33 F_q_2 = U_q * D_q * transpose(V_q);

    /////////////////////////////////////////////////// Denormalize ///////////////////////////////////////////////////

    Matrix33 T(scale,0,-scale*mean_u,
               0, scale, -scale*mean_v,
               0, 0, 1);
    Matrix33 T1(scale1,0,-scale1*mean_u1,
                0, scale1, -scale1*mean_v1,
                0, 0, 1);

    Matrix33 F = transpose(T1) * F_q_2 * T;


    /////////////////////////////////////////////// Getting the E matrix /////////////////////////////////////////////

    ///Building the K matrix:
    Matrix33 K(fx, s, cx,
               0, fy, cy,
               0, 0, 1);

    Matrix E = transpose(K) * F * K;


    /////////////////////////////////////////////// Recover R and t //////////////////////////////////////////////////




    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // TODO: Don't forget to
    //          - write your recovered 3D points into 'points_3d' (so the viewer can visualize the 3D points for you);
    //          - write the recovered relative pose into R and t (the view will be updated as seen from the 2nd camera,
    //            which can help you check if R and t are correct).
    //       You must return either 'true' or 'false' to indicate whether the triangulation was successful (so the
    //       viewer will be notified to visualize the 3D points and update the view).
    //       There are a few cases you should return 'false' instead, for example:
    //          - function not implemented yet;
    //          - input not valid (e.g., not enough points, point numbers don't match);
    //          - encountered failure in any step.
    return points_3d.size() > 0;
}

