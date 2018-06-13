/*
 * 10_LinearShooting.cpp
 * Russell Gavin Fielder, CSU Chico, Math 461 Numerical Analysis
 * Prof. Sergei Fomin
 * Version JERRY RIGGED FOR EULERS 
 *
 * This program applies the Shooting Method to a second order linear BVP
 */

/** 
 *         Instructions for Use
 *
 *  0. Ensure files "vec.h", "vec.cpp" are present
 *  1. Enter BVP Constants a, b, LBC, RBC starting at line 41
 *  2. Enter step size (h) on line 46
 *  3. Enter desired precision of output (places after decimal) on line 48
 *  4. Define the function f for y" = f(x, y, y') starting on line 195
 *  5. Compile. (For g++, the only command is "make")
 *  6. Run the output ("./linearshooting"). 
*/

#define TWO_DIM_SYSTEM_ONLY 2
#include <iostream>
#include <cmath>
#include <iomanip>
#include "vec.h"    // class vec holds n-dimensional numerical vectors of long double type
                    // also has class veclist to store a list of vec's

using namespace std;


//Function prototypes
void getNewMeshPoint(long double t, vec& y, vec& newMeshPoint, bool nonhomo);
void getNewMeshPoint_eulers(long double t, vec& y, vec& newMeshPoint, bool nonhomo);
long double f(int index, long double t, vec& y, long double addToY[], long double scalar, bool nonhomo);
long double p(long double x);
long double q(long double x);
long double r(long double x);

//BVP Constants
const int n = TWO_DIM_SYSTEM_ONLY; //Legacy constant from rungekutta.cpp, do not change
const long double a = 0.00;        
const long double b = 1.00;
const long double LBC = 0.00;      //  Left boundary condition y(a) 
const long double RBC = 0.00;      // Right boundary condition y(b)
//Step size
const long double h = 0.05;
//Program settings
const int PRECISION = 8;

//Main function
int main() {

    //Welcome message
    cout << endl;
    cout << "**********************************************\n";
    cout << "****** Linear Shooting Method results: *******\n";
    cout << "**********************************************\n\n";
    cout << "Approximating the solution of second-order BVP\n";
    cout << "with y(" << a << ") = " << LBC << ", y(" << b << ") = " << RBC << ", with h = " << h << "\n";
    cout << "and  " << a << " <= t <= " << b << " :\n\n";

    //Set precision of output
    cout << fixed << setprecision(PRECISION);
 
    //Find number of mesh points. Done this way to resolve division rounding issues.
    int numMeshPoints = 1;
    while (a + h * (numMeshPoints - 1) < b) 
        numMeshPoints++;

    // allocate memory
    veclist w1(numMeshPoints, n);
    veclist w2(numMeshPoints, n);
    long double* w = new long double [numMeshPoints];

    // Set up initial conditions -- see page 673 in text for method
    long double ICforY1[] = {LBC, 0.00};
    long double ICforY2[] = {0.00, 1.00};
    vec firstMeshPointInW1(ICforY1, n);
    vec firstMeshPointInW2(ICforY2, n);
    w1[0] = firstMeshPointInW1;
    w2[0] = firstMeshPointInW2;
   
    //Generate mesh points for w1
    int i = 0;
    while (a + i*h < b) {
        getNewMeshPoint_eulers(a + i*h, w1[i], w1[i+1], true);
        i++;
    }
    //Generate mesh points for w2
    i = 0;
    while (a + i*h < b) {
        getNewMeshPoint_eulers(a + i*h, w2[i], w2[i+1], false);
        i++;
    } 
    
    //Form the final solution approximation
    //Find error coefficient EC -- see page 673 in text for method
    long double EC = (RBC - w1[numMeshPoints-1].e[0]) / (w2[numMeshPoints-1].e[0]);
    //Define solution mesh 
    for (i = 0; i < numMeshPoints; i++) 
        w[i] = w1[i].e[0] + (EC * w2[i].e[0]);

    //Print mesh points as a table
    for (i = 0; i < numMeshPoints; i++) {   
        cout << "t = " << a + i*h << " \t>>\t";
        cout << "w" << " = " << w[i] << endl;
    }
    cout << endl;

    //Print mesh points as a wolfram-language list
    cout << "For Mathematica input:\n";
    cout << "w = {";
    for (i = 0; i < numMeshPoints - 1; i++) 
        cout << w[i] << ", ";
    cout << w[i] << "};" << endl;
    //Print t values as a wolfram-language list
    cout << "t = {";
    for (i = 0; i < (numMeshPoints - 1); i++)
        cout << (a + i*h) << ", ";
    cout << (a + i*h) << "};" << endl << endl;

    //If requested, Print mesh points of w1, w2 used in approximation
    cout << "Would you like to print the y1 and y2 approximations? y/n\n";
    char userChar;
    cin >> userChar;
    cout << endl;
    if (userChar == 'y') {
        cout << "Recall that w1 approximates {y, y'} for nonhomogenous system\n";
        cout << "Recall that w2 approximates {y, y'} for homogenous system\n";
        for (i = 0; i < numMeshPoints; i++) {
           //Print mesh points as a table
           cout << "t = " << a + i*h << " \t>>\t";
           cout << "w1" << " = ";
           w1[i].print(); 
           cout << "\t w2 = ";
           w2[i].print();
           cout << endl;
        }
        cout << endl;
    }
    
    //Deallocate memory
    w1.clear();
    w2.clear();
    delete[] w;
    firstMeshPointInW1.deallocate();
    firstMeshPointInW2.deallocate();

    return 0;
}

//uses eulers instead

void getNewMeshPoint_eulers(long double t, vec& y, vec& newMeshPoint, bool nonhomo) {
    long double k0[n]; //array of zeroes for passing a null addToY for computing k1's
    for (int i = 0; i < n; i++)
        k0[i] = 0.00;
    for (int i = 0; i < n; i++)
        newMeshPoint[i] = y[i] + (h * f(i, t, y, k0, 0.00, nonhomo));


}

/*
//This function calculates a new mesh point and stores it in a passed vec reference
void getNewMeshPoint(long double t, vec& y, vec& newMeshPoint) {
    //Temp array to hold function results
    for (int i = 0; i < n; i++) {
        newMeshPoint[i] = y[i] + (h * f[i]);
    }
}
*/

//This function calculates a new mesh point and stores it in a passed vec reference
//nonhomo = true gets the meshpoint for y1, false gets meshpoint for y2
void getNewMeshPoint(long double t, vec& y, vec& newMeshPoint, bool nonhomo) {
    long double k0[n]; //array of zeroes for passing a null addToY for computing k1's
    for (int i = 0; i < n; i++)
        k0[i] = 0.00;
    long double k1[n];
    long double k2[n];
    long double k3[n];
    long double k4[n];
    //Note: all the k1's must be computed before any of the k2's can be computed
    for (int i = 0; i < y.dim; i++) 
        k1[i] = h * f(i, t, y, k0, 0.00, nonhomo);
    //All the k1's computed, now we can compute the k2's
    for (int i = 0; i < n; i++)
        k2[i] = h * f(i, t + (h/2), y, k1, 0.50, nonhomo);
    //All the k2's computed, now we can compute the k3's
    for (int i = 0; i < n; i++)
        k3[i] = h * f(i, t + (h/2), y, k2, 0.50, nonhomo);
    //All the k3's computed, now we can compute the k4's
    for (int i = 0; i < n; i++)
        k4[i] = h * f(i, t + h, y, k3, 1.00, nonhomo);
    //Now lets construct the meshpoint
    for (int i = 0; i < n; i++)
        newMeshPoint[i] = y[i] + ((k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.00);

}

/**
 * This function evaluates the vector function f
 *
 * @param index           which component of f to evaluate
 * @param x               the x value to evaluate the function component at
 * @param u_unmodded      a vec reference storing the previous meshpoint
 * @param addToU[]        an array of long doubles to add to the components of 
 *                            the previous meshpoint. This is k1[i], k2[i], etc.
 * @param scalar          a scalar to multiply each addToU value before adding
 * @param nonhomo         whether to add the nonhomogeneous portion r(x) (for y1 solution)
 *
 * @return                the value of the requested component of f
*/
long double f(int index, long double x, vec& u_unmodded, long double addToU[], long double scalar, bool nonhomo) {
    long double result = 0.00;
    long double u[n];
    for (int i = 0; i < n; i++)
        u[i] = u_unmodded[i] + (scalar * addToU[i]);
    switch(index) {
        case 0:
            // First state variable derivative -- do not change (f definition is below)
            result = u[1];
            break;
        case 1:
            // Second state variable derivative -- do not change (f definition is below)
            result = p(x)*u[1] + q(x)*u[0];
            if (nonhomo) result += r(x); 
            break;
    }
    return result;
}


//************************************************************************
//************************ Define f(x, y, y') Here ***********************
//************************************************************************

//this function defines p(x), the coefficient on y'  in y'' = p(x)y' + q(x)y + r(x)
long double p(long double x) {

    long double P_x = 5.00;

    return P_x;
}
//this function defines q(x), the coefficient on y  in y'' = p(x)y' + q(x)y + r(x)
long double q(long double x) {

    long double Q_x = 0.00;

    return Q_x;
}
//this function defines r(x), the nonhomogeneous term in y'' = p(x)y' + q(x)y + r(x)
long double r(long double x) {

    long double R_x = 3*x*x - sin(5*x);
    
    return R_x;
}

//************************************************************************
//**************************** End Define f ******************************
//************************************************************************









