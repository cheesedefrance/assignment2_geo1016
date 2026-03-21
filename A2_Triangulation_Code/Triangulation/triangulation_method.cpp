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

// helper function to triangulate single point pair given two camera projections
Vector3D triangulate_points(
        const Vector2D& p0, // observed points in image 1
        const Vector2D& p1, // observed points in image 2
        const Matrix34& M0, // camera matrix 1 (K[I|0])
        const Matrix34& M1) // camera matrix 2 (K[R|t])
{
    double x0 = p0[0], y0 = p0[1];
    double x1 = p1[0], y1 = p1[1];


    Matrix A(4, 4);
    A.set_row(0, x0*M0.get_row(2) - M0.get_row(0));
    A.set_row(1, y0*M0.get_row(2) - M0.get_row(1));
    A.set_row(2, x1*M1.get_row(2) - M1.get_row(0));
    A.set_row(3, y1*M1.get_row(2) - M1.get_row(1));

    // solving for AP = 0 via SVD with P is last column of V
    Matrix U4(4, 4, 0.0), S4(4, 4, 0.0), V4(4, 4, 0.0);
    svd_decompose(A, U4, S4, V4);

    Vector hp3 = V4.get_column(V4.cols() - 1); // homogeneous 3D point

    // divide X,Y,Z by W to get unhomogeneous coordinates
    return Vector3D(hp3[0] / hp3[3],
                    hp3[1] / hp3[3],
                    hp3[2] / hp3[3]);
}
//helper function to calculate centroid
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
//helper function to calculate scale factor s
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


    int n = points_0.size(); //variable used in for loops
    // calculating centroid for point sets 0 and 1
    Vector2D c0 = get_centroid(points_0);
    Vector2D c1 = get_centroid(points_1);

    std::cout << "centroid points_0 -- " << c0 << "\ncentroid points_1 -- " << c1 << "\n" << std::endl;

    // calculating scaling factor for psets 0 and 1
    double s0 = get_s(points_0, c0);
    double s1 = get_s(points_1, c1);

    // making normalization matrix T for psets 0 and 1
    Matrix33 T0(s0, 0, - s0 * c0[0],
        0, s0, - s0 * c0[1],
        0, 0, 1);

    Matrix33 T1(s1, 0, - s1 * c1[0],
        0, s1,  - s1 * c1[1],
        0, 0, 1);


    std::cout << "Matrix T0 -- " << T0 << "\nMatrix T1 -- " << T1 << "\n" << std::endl;


    // applying normalization to psets 0 and 1 by tnp = T * np
    std::vector<Vector3D> np0(points_0.size());
    std::vector<Vector3D> np1(points_1.size());


    for (int i=0; i<n; i++) { //compute normalized points
        np0[i] = T0*points_0[i].homogeneous();
        np1[i] = T1*points_1[i].homogeneous();
    }
    // initializing matrix W
    Matrix W(n, 9);

    for (int i = 0; i < n; i++) {
        double u0  = np0[i].x();
        double v0  = np0[i].y();
        double u1 = np1[i].x();
        double v1 = np1[i].y();

        W.set_row(i, {u0*u1, v0*u1, u1, u0*v1, v0*v1, v1, u0, v0, 1}); //populate matrix W
    }

    Matrix U(n, n, 0.0), S(n, 9, 0.0), V(9, 9, 0.0);
    svd_decompose(W,U,S,V);

    // calculating the fundamental matrix F (slide 46 and notes p.8). this one is still normalized, thus fh or f hat name
    Vector fh = V.get_column(V.cols() - 1);
    Matrix FH(3, 3, fh.data());

    // enforcing rank(F) = 2
    Matrix U2(3, 3, 0.0), S2(3, 3, 0.0), V2(3, 3, 0.0);
    svd_decompose(FH,U2,S2,V2);

    S2(2, 2) = 0.0;
    Matrix33 FQ = U2 * S2 * transpose(V2); //reconstruct FQ, fundamental matrix with enforced rank=2 but still normalized
    Matrix33 F = transpose(T1)*FQ*T0; //recover unnormalized F

    std::cout << "fundamental matrix FH\n" <<  FH << "\nfundamental matrix FQ\n" << FQ << "\nfundamental matrix F\n" << F  << std::endl;

    // computing essential matrix E
    Matrix33 K(fx, s, cx, 0, fy, cy, 0, 0, 1); //construct matrix K from given camera parameters
    Matrix33 E = transpose(K) * F * K; // denormalized F

    // recovering R, t from SVD(E)
    Matrix U3(3, 3, 0.0), S3(3, 3, 0.0), V3(3, 3, 0.0);
    svd_decompose(E, U3, S3, V3);

    Matrix33 WE(0, -1, 0, 1, 0, 0, 0, 0, 1);

    Vector t1 = U3.get_column(2);
    Vector t2 = - U3.get_column(2);

    Matrix33 R1 = determinant(U3 * WE * transpose(V3)) * U3 * WE * transpose(V3); //enforce determinant positive 1
    Matrix33 R2 = determinant(U3 * transpose(WE) * transpose(V3)) * U3 * transpose(WE) * transpose(V3);

    std::cout << "t1 -- " << t1 << "\n" << "t2 -- " << t2 << std::endl;
    std::cout << "\nR1 -- " << R1 << "\n" << "R2 -- " <<  R2 << "\n" << std::endl;

    std::cout << "(2) fundamental matrix F, essential matrix E and R, t computed.\n" << std::endl;


    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // P0 camera #1 projection matrix formula K[I|0]
    Matrix34 M0(fx,  s, cx, 0,
                 0, fy, cy, 0,
                 0,  0,  1, 0);

    std::vector<Matrix33> Rpos = { R1, R1, R2, R2 };
    std::vector<Vector3D> tpos = { t1, t2, t1, t2 };

    int bestcount = -1;
    int bestindex =  0;

    for (int c = 0; c < 4; c++) {
        Matrix33 Rc = Rpos[c];
        Vector3D tc = tpos[c];

        // // M1 camera 2 projection matrix, from K[R|t]

        Matrix33 KR = K * Rc;
        Vector3D Kt = K * tc;
        Matrix34 M1(
            KR(0,0), KR(0,1), KR(0,2), Kt[0],
            KR(1,0), KR(1,1), KR(1,2), Kt[1],
            KR(2,0), KR(2,1), KR(2,2), Kt[2]
        );

        int count = 0;
        for (int i = 0; i < n; i++) {
            Vector3D P = triangulate_points(points_0[i], points_1[i], M0, M1);

            // point must have positive depth in camera #1 (Z > 0)
            bool infront0 = P.z() > 0;

            // transform point into camera 2's frame and depth in camera #2 (Z > 0)
            Vector3D  transformedP = Rc * P + tc;
            bool infront1 = transformedP.z() > 0;

            if (infront0 && infront1) //if point has z>0 for both camera views, add one to count
                count++;
        }

        std::cout << "hypothesis " << c << " -- " << count << " points in front.\n";

        // bestcount = -1 so this check will always run
        if (count > bestcount) {
            bestcount = count;
            bestindex = c;
        }
    }

    if (bestcount==-1) { //in case no points lie in front of cameras, abort program
        std::cout<< "No combination of camera views is correct" <<std::endl;
        return false;
    }

    // pick the best relative pose (R, t)
    R = Rpos[bestindex];
    t = tpos[bestindex];

    std::cout << "\nbest relative pose -- " << bestindex
              << "\nwith " << bestcount << " / " << n << " points in front\n" << std::endl;

    std::cout << "final R -- " << R << "\nfinal t -- " << t << "\n" << std::endl;

    std::cout << "(3) best R, t selected.\n" << std::endl;

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

    // winning relative pose M_final
    Matrix33 KR = K * R;
    Vector3D Kt = K * t;
    Matrix34 M_final(
    KR(0,0), KR(0,1), KR(0,2), Kt[0],
    KR(1,0), KR(1,1), KR(1,2), Kt[1],
    KR(2,0), KR(2,1), KR(2,2), Kt[2]
    );


    for (int i = 0; i < n; i++) {
        Vector3D P = triangulate_points(points_0[i], points_1[i], M0, M_final);
        points_3d.push_back(P);
    }

    std::cout << "(4) 3D points reconstructed.\n\n" << std::flush;


    //checking reprojection error?
    double reprojection_error = 0.0;
    for (int i = 0; i < n; i++) {
        // project into image 1
        Vector3D p1 = K * points_3d[i];
        double u1 = p1.x()/p1.z(), v1 = p1.y()/p1.z();

        // project into image 2
        Vector3D p2 = K * (R * points_3d[i] + t);
        double u2 = p2.x()/p2.z(), v2 = p2.y()/p2.z();

        double e1 = (u1-points_0[i].x())*(u1-points_0[i].x())
                  + (v1-points_0[i].y())*(v1-points_0[i].y());
        double e2 = (u2-points_1[i].x())*(u2-points_1[i].x())
                  + (v2-points_1[i].y())*(v2-points_1[i].y());

        reprojection_error += sqrt(e1) + sqrt(e2);
    }
    reprojection_error /= (2*n);
    std::cout << "Mean reprojection error: " << reprojection_error << " px" << std::endl;
    return points_3d.size() > 0;
}