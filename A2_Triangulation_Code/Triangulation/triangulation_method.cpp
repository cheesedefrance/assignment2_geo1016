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
        const Vector2D& p0, // observed point in image 1
        const Vector2D& p1, // observed point in image 2
        const Matrix34& P0, // camera matrix 1 (K[I|0])
        const Matrix34& P1) // camera matrix 2 (K[R|t])
{
    double x0 = p0[0], y0 = p0[1];
    double x1 = p1[0], y1 = p1[1];

    // building matrix A (each view contributes two rows)
    Matrix A(4, 4);

    // view #1 rows;
    // x*(m3^T) - m1^T
    // y*(m3^T) - m2^T
    A.set_row(0, { x0*P0(2,0) - P0(0,0),
                   x0*P0(2,1) - P0(0,1),
                   x0*P0(2,2) - P0(0,2),
                   x0*P0(2,3) - P0(0,3) });

    A.set_row(1, { y0*P0(2,0) - P0(1,0),
                   y0*P0(2,1) - P0(1,1),
                   y0*P0(2,2) - P0(1,2),
                   y0*P0(2,3) - P0(1,3) });

    // view #2 rows;
    // x'*(m3'^T) - m1'^T
    // y'*(m3'^T) - m2'^T
    A.set_row(2, { x1*P1(2,0) - P1(0,0),
                   x1*P1(2,1) - P1(0,1),
                   x1*P1(2,2) - P1(0,2),
                   x1*P1(2,3) - P1(0,3) });

    A.set_row(3, { y1*P1(2,0) - P1(1,0),
                   y1*P1(2,1) - P1(1,1),
                   y1*P1(2,2) - P1(1,2),
                   y1*P1(2,3) - P1(1,3) });

    // solving for AP = 0 via SVD with P is last column of V
    Matrix U4(4, 4, 0.0), S4(4, 4, 0.0), V4(4, 4, 0.0);
    svd_decompose(A, U4, S4, V4);

    Vector hp3 = V4.get_column(V4.cols() - 1); // homogeneous 3D point

    // divide X,Y,Z by W to get unhomogeneous coordinates
    return Vector3D(hp3[0] / hp3[3],
                    hp3[1] / hp3[3],
                    hp3[2] / hp3[3]);
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

    Vector t1 = U3.get_column(U3.cols() - 1);
    Vector t2 = - U3.get_column(U3.cols() - 1);

    Matrix33 R1 = determinant(U3 * WE * transpose(V3)) * U3 * WE * transpose(V3);
    Matrix33 R2 = determinant(U3 * transpose(WE) * transpose(V3)) * U3 * transpose(WE) * transpose(V3);

    std::cout << "t1 -- " << t1 << "\n" << "t2 -- " << t2 << std::endl;
    std::cout << "\nR1 -- " << R1 << "\n" << "R2 -- " <<  R2 << "\n" << std::endl;

    std::cout << "(2) fundamental matrix F, essential matrix E and R, t computed.\n" << std::endl;


    // TODO: Reconstruct 3D points. The main task is
    //      - triangulate a pair of image points (i.e., compute the 3D coordinates for each corresponding point pair)

    // P0 camera #1 projection matrix formula K[I|0]
    Matrix34 P0(fx,  s, cx, 0,
                 0, fy, cy, 0,
                 0,  0,  1, 0);

    std::vector<Matrix33> Rpos = { R1, R1, R2, R2 };
    std::vector<Vector3D> tpos = { t1, t2, t1, t2 };

    int bestcount = -1;
    int bestindex =  0;

    for (int c = 0; c < 4; c++) {
        Matrix33 Rc = Rpos[c];
        Vector3D tc = tpos[c];

        // P1 camera #2 projection matrix formula K[R|t]
        Matrix34 P1(
            fx * Rc(0,0) + cx * Rc(2,0), fx * Rc(0,1) + cx * Rc(2,1), fx * Rc(0,2) + cx * Rc(2,2), fx * tc[0] + cx * tc[2],
            fy * Rc(1,0) + cy * Rc(2,0),  fy * Rc(1,1) + cy * Rc(2,1),  fy * Rc(1,2) + cy * Rc(2,2),  fy * tc[1] + cy * tc[2],
            Rc(2,0), Rc(2,1), Rc(2,2), tc[2]);

        int count = 0;
        for (int i = 0; i < n; i++) {
            Vector3D p3 = triangulate_points(points_0[i], points_1[i], P0, P1);

            // point must have positive depth in camera #1 (Z > 0)
            bool infront0 = p3[2] > 0;

            // transform point into camera 2's frame and depth in camera #2 (Z > 0)
            Vector3D  transform = Rc * p3 + tc;
            bool infront1 = transform[2] > 0;

            if (infront0 && infront1)
                count++;
        }

        std::cout << "hypothesis " << c << " -- " << count << " points in front.\n";

        // bestcount = -1 so this check will always run
        if (count > bestcount) {
            bestcount = count;
            bestindex = c;
        }
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

    // winning relative pose P2
    Matrix34 P2(
    fx * R(0,0) + cx * R(2,0), fx * R(0,1) + cx * R(2,1), fx*R(0,2) + cx * R(2,2), fx * t[0] + cx * t[2],
    fy * R(1,0) + cy * R(2,0), fy * R(1,1) + cy * R(2,1), fy * R(1,2) + cy * R(2,2), fy * t[1] + cy * t[2],
    R(2,0), R(2,1), R(2,2), t[2]);

    points_3d.clear();
    points_3d.reserve(n);

    for (int i = 0; i < n; i++) {
        Vector3D p3 = triangulate_points(points_0[i], points_1[i], P0, P2);
        points_3d.push_back(p3);
    }

    std::cout << "(4) 3D points reconstructed.\n\n" << std::flush;

    return points_3d.size() > 0;
}