//
// Created by jan on 4/5/19.
//

#ifndef MAGNETICSPHERES_DIPOLE_H
#define MAGNETICSPHERES_DIPOLE_H

#include <vector>

class Dipole {
private:
    std::vector<double> r, m;
public:
    Dipole();
    void SetR(double x, double y, double z);
    void SetM(double phi, double theta, double m = 1);
    std::vector<double> GetR();
    std::vector<double> GetM();

    // mathematical properties
    double DistanceTo(Dipole v);
    std::vector<double> VectorTo(Dipole v);

    //output
    void Print();

};


#endif //MAGNETICSPHERES_DIPOLE_H
