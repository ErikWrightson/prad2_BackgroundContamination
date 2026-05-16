/**
 * Contains Utility functions that can be used independently.
 * 
 * @author Erik Wrightson <wrightso@jlab.org>
 * @version 05.14.2026
 * @creation 05.14.2026
 */
#include "Utils.h"

 Utils::LineOfBestFit Utils::FitLine(vector<Float_t> x, vector<Float_t> y, vector<Float_t> z){
    Utils::LineOfBestFit line;

    int n = x.size();

    // TPrincipal does PCA
    TPrincipal pca(3, "");

    Double_t point[3];

    for (int i = 0; i < n; i++) {
        point[0] = x[i];
        point[1] = y[i];
        point[2] = z[i];

        pca.AddRow(point);
    }

    pca.MakePrincipals();

    // Mean point: point on the fitted line
    const TVectorD *mean = pca.GetMeanValues();

    line.x0 = (*mean)[0];
    line.y0 = (*mean)[1];
    line.z0 = (*mean)[2];

    // Eigenvectors: first principal component gives best-fit direction
    const TMatrixD *eig = pca.GetEigenVectors();

    line.vx = (*eig)(0,0);
    line.vy = (*eig)(1,0);
    line.vz = (*eig)(2,0);

    return line;
}

/**
 * Returns the point of the closest approach from the z axis for given the line of best fit handed in.
 *
 * @param line - the line of best fit to find its approach to the z-axis.
 *
 * @return - the point of closest approach to the z axis.
 */
Utils::Point Utils::ClosestApproachToZAxis(Utils::LineOfBestFit line){
    Float_t t = -1.0 * ((line.x0*line.vx)+(line.y0*line.vy))/((line.vx*line.vx) + (line.vy*line.vy));

    Utils::Point p;
    p.x = line.x0 + t*line.vx;
    p.y = line.y0 + t*line.vy;
    p.z = line.z0 + t*line.vz;

    return p;
}