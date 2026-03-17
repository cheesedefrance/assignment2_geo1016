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

#include <array>

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

    // std::cout << "\nTODO: implement the 'triangulation()' function in the file 'Triangulation/triangulation_method.cpp'\n\n";
    //
    // std::cout << "[Liangliang]:\n"
    //              "\tSimilar to the first assignment, basic linear algebra data structures and functions are provided in\n"
    //              "\tthe following files:\n"
    //              "\t    - Triangulation/matrix.h: handles matrices of arbitrary dimensions and related functions.\n"
    //              "\t    - Triangulation/vector.h: manages vectors of arbitrary sizes and related functions.\n"
    //              "\t    - Triangulation/matrix_algo.h: contains functions for determinant, inverse, SVD, linear least-squares...\n"
    //              "\tFor more details about these data structures and a complete list of related functions, please\n"
    //              "\trefer to the header files mentioned above.\n\n"
    //              "\tIf you choose to implement the non-linear method for triangulation (optional task). Please\n"
    //              "\trefer to 'Tutorial_NonlinearLeastSquares/main.cpp' for an example and some explanations.\n\n"
    //              "\tFor your final submission, adhere to the following guidelines:\n"
    //              "\t    - submit ONLY the 'Triangulation/triangulation_method.cpp' file.\n"
    //              "\t    - remove ALL unrelated test code, debugging code, and comments.\n"
    //              "\t    - ensure that your code compiles and can reproduce your results WITHOUT ANY modification.\n\n" << std::flush;
    //
    // /// Below are a few examples showing some useful data structures and APIs.
    //
    // /// define a 2D vector/point
    // Vector2D b(1.1, 2.2);
    //
    // /// define a 3D vector/point
    // Vector3D a(1.1, 2.2, 3.3);
    //
    // /// get the Cartesian coordinates of a (a is treated as Homogeneous coordinates)
    // Vector2D p = a.cartesian();
    //
    // /// get the Homogeneous coordinates of p
    // Vector3D q = p.homogeneous();
    //
    // /// define a 3 by 3 matrix (and all elements initialized to 0.0)
    // Matrix33 A;
    //
    // /// define and initialize a 3 by 3 matrix
    // Matrix33 T(1.1, 2.2, 3.3,
    //            0, 2.2, 3.3,
    //            0, 0, 1);
    //
    // /// define and initialize a 3 by 4 matrix
    // Matrix34 M(1.1, 2.2, 3.3, 0,
    //            0, 2.2, 3.3, 1,
    //            0, 0, 1, 1);
    //
    // /// set first row by a vector
    // M.set_row(0, Vector4D(1.1, 2.2, 3.3, 4.4));
    //
    // /// set second column by a vector
    // M.set_column(1, Vector3D(5.5, 5.5, 5.5));
    //
    // /// define a 15 by 9 matrix (and all elements initialized to 0.0)
    // Matrix W(15, 9, 0.0);
    // /// set the first row by a 9-dimensional vector
    // W.set_row(0, {0, 1, 2, 3, 4, 5, 6, 7, 8}); // {....} is equivalent to a std::vector<double>
    //
    // /// get the number of rows.
    // int num_rows = W.rows();
    //
    // /// get the number of columns.
    // int num_cols = W.cols();
    //
    // /// get the the element at row 1 and column 2
    // double value = W(1, 2);
    //
    // /// get the last column of a matrix
    // Vector last_column = W.get_column(W.cols() - 1);
    //
    // /// define a 3 by 3 identity matrix
    // Matrix33 I = Matrix::identity(3, 3, 1.0);
    //
    // /// matrix-vector product
    // Vector3D v = M * Vector4D(1, 2, 3, 4); // M is 3 by 4
    //
    // ///For more functions of Matrix and Vector, please refer to 'matrix.h' and 'vector.h'

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function)
    // checking valid input.
    if (points_0.size() < 8 || points_1.size() < 8){
        std::cout << "invalid input. point size smaller than eight." << std::endl;
        return false;
    }
    if (points_0.size() != points_1.size()){
        std::cout << "invalid input. different point sizes." << std::endl;
        return false;
    }

    std::cout << "\n(1) input checked.\n" << std::endl;


    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    // calculating centroid for point sets 0 and 1
    double mx0 = 0, my0 = 0;
    for (auto& p: points_0) {
        mx0 += p[0];
        my0 += p[1];

        // std::cout << mx << my << std::endl;
    }

    Vector2D c0(mx0 / points_0.size(), my0 / points_0.size());

    double mx1 = 0, my1 = 0;
    for (auto& p: points_1) {
        mx1 += p[0];
        my1 += p[1];

        // std::cout << mx << my << std::endl;
    }

    Vector2D c1(mx1 / points_1.size(), my1 / points_1.size());

    std::cout << "centroid points_0 -- " << c0 << "\ncentroid points_1 -- " << c1 << "\n" << std::endl;

    // calculating average distance for psets 0 and 1
    double dsum0 = 0;
    for (auto& p : points_0) {
        double dx = p[0] - c0[0];
        double dy = p[1] - c0[1];
        dsum0 += std::sqrt(pow(dx, 2) + pow(dy, 2));
    }

    double d0 = dsum0 / points_0.size();

    double dsum1 = 0;
    for (auto& p : points_1) {
        double dx = p[0] - c1[0];
        double dy = p[1] - c1[1];
        dsum1 += std::sqrt(pow(dx, 2) + pow(dy, 2));
    }

    double d1 = dsum1 / points_1.size();

    std::cout << "d0 -- " << d0 << "\nd1 -- " << d1 << "\n" << std::endl;

    // calculating scaling factor for psets 0 and 1
    double s0 = sqrt(2) / d0;
    double s1 = sqrt(2) / d1;

    // making normalization matrix T for psets 0 and 1
    Matrix33 T0(s0, 0, - s0 * c0[0],
        0, s0, - s0 * c0[1],
        0, 0, 1);

    Matrix33 T1(s1, 0, - s1 * c1[0],
        0, s1,  - s1 * c1[1],
        0, 0, 1);


    std::cout << "Matrix T0 -- " << T0 << "\nMatrix T1 -- " << T1 << "\n" << std::endl;

    // making psets 0 and 1 into homogeneous coordinates
    std::vector<Vector3D> hp0;
    hp0.resize(points_0.size());
    for (size_t i = 0; i < points_0.size(); i++) {
        hp0[i] = Vector3D(points_0[i][0], points_0[i][1], 1);
    }

    std::vector<Vector3D> hp1;
    hp1.resize(points_1.size());
    for (size_t i = 0; i < points_1.size(); i++) {
        hp1[i] = Vector3D(points_1[i][0], points_1[i][1], 1);
    }

    // applying normalization to psets 0 and 1 by tnp = T * np
    std::vector<Vector3D> np0(hp0.size());
    for (size_t i = 0; i < hp0.size(); i++){
        np0[i] = T0 * hp0[i];
    }

    std::vector<Vector3D> np1(hp1.size());
    for (size_t i = 0; i < hp1.size(); i++){
        np1[i] = T1 * hp1[i];
    }

    // making matrix W
    int n = np0.size();
    Matrix W(n, 9);

    for (int i = 0; i < n; i++) {
        double u0  = np0[i][0];
        double v0  = np0[i][1];
        double u1 = np1[i][0];
        double v1 = np1[i][1];

        W.set_row(i, {u0*u1, v0*u1, u1, u0*v1, v0*v1, v1, u0, v0, 1});
    }

    Matrix U(n, n, 0.0);
    Matrix S(n, 9, 0.0);
    Matrix V(9, 9, 0.0);
    svd_decompose(W,U,S,V);

    // std::cout << "USV check." << std::endl;
    // // check 1: U is orthogonal, so U * U^T must be identity
    // std::cout << "U*U^T: \n" << U * transpose(U) << std::endl;
    // // check 2: V is orthogonal, so V * V^T must be identity
    // std::cout << "V*V^T: \n" << V * transpose(V) << std::endl;
    // // check 3: S must be a diagonal matrix
    // std::cout << "S: \n" << S << std::endl;

    // calculating the fundamental matrix F (slide 46 and notes p.8)
    Vector fh = V.get_column(V.cols() - 1);
    Matrix FH(3, 3);
    FH.set_row(0, {fh[0], fh[1], fh[2]});
    FH.set_row(1, {fh[3], fh[4], fh[5]});
    FH.set_row(2, {fh[6], fh[7], fh[8]});

    // enforcing rank(F) = 2
    Matrix U2(3, 3, 0.0);
    Matrix S2(3, 3, 0.0);
    Matrix V2(3, 3, 0.0);
    svd_decompose(FH,U2,S2,V2);

    S2(2, 2) = 0.0;
    Matrix33 FQ = U2 * S2 * transpose(V2);
    Matrix33 F = transpose(T1)*FQ*T0;

    std::cout << U2 << S2 << V2 << std::endl;
    std::cout << "fundamental matrix FH\n" <<  FH << "\nfundamental matrix FQ\n" << FQ << "\nfundamental matrix F\n" << F  << std::endl;

    // computing essential matrix E
    Matrix33 K(fx, s, cx, 0, fy, cy, 0, 0, 1);
    Matrix33 E = transpose(K) * F * K; // denormalized F

    // recovering R, t from SVD(E)
    Matrix U3(3, 3, 0.0);
    Matrix S3(3, 3, 0.0);
    Matrix V3(3, 3, 0.0);
    svd_decompose(E, U3, S3, V3);

    Matrix33 WE(0, -1, 0, 1, 0, 0, 0, 0, 1);
    Matrix33 ZE(0, 1, 0, -1, 0, 0, 0, 0, 0);

    std::cout << U3 << S3 << V3 << std::endl;

    Vector t1 = U3.get_column(U3.cols() - 1);
    Vector t2 = - U3.get_column(U3.cols() - 1);

    Matrix33 R1 = U3 * WE * transpose(V3);
    Matrix33 R2 = U3 * transpose(WE) * transpose(V3);

    if (determinant(R1) < 0){
        R1 = - R1;
    }
    if (determinant(R2) < 0){
        R2 = - R2;
    }

    std::cout << "(2) fundamental matrix F, essential matrix E and R, t computed." << std::endl;


    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // four possible combinations for R, t
    std::vector<Matrix33> Rpos = {R1, R1, R2, R2};
    std::vector<Vector3D> tpos = {t1, t2, t1, t2};

    int bestcount = -1;
    int bestindex = 0;

    // for (int c = 0; c < 4; c++){
    //     Matrix33 R = Rpos[c];
    //     Vector3D t = tpos[c];
    //
    //     // camera #1
    //     Matrix34 P0(fx, s, cx, 0, 0, fy, cy, 0, 0, 0, 1, 0);
    //
    //     // camera #2
    //     Matrix34 P1(
    //     fx*(R(0,0)) + cx*R(2,0), fx*(R(0,1)) + cx*R(2,1), fx*(R(0,2)) + cx*R(2,2), fx*t[0] + cx*t[2],
    //     fy*(R(1,0)) + cy*R(2,0), fy*(R(1,1)) + cy*R(2,1), fy*(R(1,2)) + cy*R(2,2), fy*t[1] + cy*t[2],
    //     R(2,0), R(2,1), R(2,2), t[2]);
    //
    //     int count = 0;
    //     for (int i = 0; i < n; i++){
    //         Vector3D P = triangulation(points_0[i], points_1[i], P0, P1);
    //
    //         // check in front of camera #1
    //         bool infront0 = P[2] > 0;
    //
    //         // check in front of camera #2
    //         Vector3D Q = R * P + t;
    //         bool infront1 = Q[2] > 0;
    //
    //         if (infront0 && infront1){
    //             count++;
    //         }
    //
    //         if (count > bestcount)
    //         {
    //             bestcount = count;
    //             bestindex = c;
    //         }
    //     }
    // }
    //
    // // correct R, t
    // R = Rpos[bestindex];
    // t = tpos[bestindex];

    std::cout << "\n(3) text.\n" << std::endl;


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

    std::cout << "\n(4) text.\n" << std::endl;

    return points_3d.size() > 0;
}