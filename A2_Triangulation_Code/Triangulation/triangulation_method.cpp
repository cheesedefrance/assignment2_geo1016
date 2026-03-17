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

Vector2D get_centroid(const std::vector<Vector2D> &points) {
    Vector2D centroid{0.0,0.0};
    for (int i = 0; i < points.size(); i++) {
        centroid.x() += points[i].x();
        centroid.y() += points[i].y();
    }

    centroid.x() /= points.size();
    centroid.y() /= points.size();

    return centroid;
}


double get_s(const std::vector<Vector2D> &points, Vector2D cent){
    double mean_distance = 0.0;

    for (int i=0; i< points.size(); i++) {
        double xdiff = points[i].x() - cent.x();
        double ydiff = points[i].y() - cent.y();
        mean_distance += sqrt(xdiff*xdiff + ydiff*ydiff);
    }

    mean_distance /= points.size();

    double s = sqrt(2.0)/mean_distance;



    return s;

}

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
    

    // TODO: delete all above example code in your final submission

    //--------------------------------------------------------------------------------------------------------------
    // implementation starts ...

    // TODO: check if the input is valid (always good because you never known how others will call your function).

    if (points_0.size() != points_1.size()) {
        std::cout << "invalid input, inputs have different size" << std::endl;
        return false;
    } else if (points_0.size() <8 || points_1.size() < 8) {
        std::cout << "invalid input, less than eight inputs" << std::endl;
        return false;
    }

    // TODO: Estimate relative pose of two views. This can be subdivided into
    //      - estimate the fundamental matrix F;
    //      - compute the essential matrix E;
    //      - recover rotation R and t.

    int n = points_0.size();

    //Matrix33 T0(get_T_matrix(points_0));
    //Matrix33 T1(get_T_matrix(points_1));


    Vector2D centroid0 = get_centroid(points_0);
    Vector2D centroid1 = get_centroid(points_1);

    double s0 = get_s(points_0, centroid0);
    double s1 = get_s(points_1, centroid1);



    std::cout << s0 << std::endl;
    std::cout << s1 << std::endl;
    std::cout << centroid0 << std::endl;
    std::cout << centroid1 << std::endl;

    std::cout << "got centroids and s" << std::endl;
    std::vector<Vector3D> normalized0;
    normalized0.resize(points_0.size());
    std::vector<Vector3D> normalized1;
    normalized1.resize(points_1.size());

    Matrix33 T0(s0, 0, -s0*centroid0.x(), 0, s0, -s0*centroid0.y(), 0, 0, 1);
    Matrix33 T1(s1, 0, -s1*centroid1.x(), 0, s1, -s1*centroid1.y(), 0, 0, 1);

    std::cout << "reserved normalized" << std::endl;

std::cout << T0 << std::endl;
std::cout << points_0[0].homogeneous() << std::endl;

    std::cout << T0*points_0[0].homogeneous() << std::endl;


    for (int i=0; i<n; i++) {
        normalized0[i] = T0*points_0[i].homogeneous();
        normalized1[i] = T1*points_1[i].homogeneous();
    }

    std::cout << "ran normalization" << std::endl;

    ///std::vector<Vector2D> normpoints0 = T0 points_0;
    Matrix W(n, 9, 0.0);

    // i++ in c++ means incrementing for each element
    for (int i = 0; i<n; i++) {
        double u = normalized0[i].x();
        double uprime = normalized1[i].x();
        double v = normalized0[i].y();
        double vprime = normalized1[i].y();

        W.set_row(i, {u*uprime, v*uprime, uprime, u*vprime, v*vprime, vprime, u, v, 1}); // using the equations from slide 25 lecture 3
    }

    std::cout << W << std::endl;

    Matrix U(n, n, 0.0);
    Matrix S(n, 9, 0.0);
    Matrix V(9, 9, 0.0);

    // calling the function to get the V (last column) which contains the smallest eigenvalues
    // to minimize the error
    svd_decompose(W,U,S,V);


    Vector fhat = V.get_column(V.cols() - 1);

    Matrix Fhat(3,3,fhat.data());

    std::cout << Fhat << std::endl;

    Matrix U2(3, 3, 0.0);
    Matrix S2(3, 3, 0.0);
    Matrix V2(3, 3, 0.0);

    svd_decompose(Fhat, U2, S2, V2);
    std::cout << "u2" << U2 << "s2" << S2 << "v2" << V2 << std::endl;
    S2(2,2) = 0;
    std::cout << "new s2" <<S2 << std::endl;
    Matrix Fq = U2*S2*transpose(V2);


    //get fundamental matrix F
    Matrix33 F = transpose(T1)*Fq*T0;
    std::cout << "F:" << F << std::endl;

    Matrix33 K(fx, s, cx, 0, fy, cy, 0, 0, 1);

    //get essential matrix E
    Matrix33 E = transpose(K)*F*K;

    std::cout << "E:" << E << std::endl;

    //get rotation and translation
    Matrix U3(3, 3, 0.0);
    Matrix S3(3, 3, 0.0);
    Matrix V3(3, 3, 0.0);

    svd_decompose(E, U3, S3, V3);
    std::cout << "U3" << U3 << "s3" << S3 << "v3" << V3 << std::endl;

    Matrix33 W2(0, -1, 0, 1, 0, 0, 0, 0, -1);

    Matrix R1 = U3*W2*transpose(V3);
    Matrix R2 = U3*transpose(W2)*transpose(V3);

    Vector3D t1= U3.get_column(2);
    Vector3D t2= -U3.get_column(2);


    std::cout<<" R1 "<< R1 << " R2 " << R2 << " t1 " << t1 << " t2 " << t2 << std::endl;
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