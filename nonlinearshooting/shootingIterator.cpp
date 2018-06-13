/*
 * shootingIterator.cpp
 * Russell Gavin Fielder, CSU Chico, Math 461 Numerical Analysis
 * Prof. Sergei Fomin
 * Version 4/30/2017
 *
 * This file is the implementation for class shootingIterator, which handles
 * the nonlinear shooting method for second-degree boundary value problems
 */

#include "shootingIterator.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;

//Pointers to the functions f and g which are defined by the user in main file
long double (*f)(long double x, long double y1, long double y2) = 0;
long double (*g)(long double x, long double z1, long double z2, 
                         long double* y_approx, long double* yprime_approx,
                         long double* y_ful, long double* yp_ful) = 0;
//Static member variables
long double * shootingIterator::x = 0;
long double * shootingIterator::y_approximation = 0;
long double * shootingIterator::yprime_approximation = 0;
long double * shootingIterator::z_approximation = 0;
long double * shootingIterator::zprime_approximation = 0;
long double * shootingIterator::y_full = 0;
long double * shootingIterator::yprime_full = 0;
int shootingIterator::selectedMeshPointFunction = 0;
int shootingIterator::selectedTValueIterator = 0;


//Constructor
// simple form: shootingIterator(a, b, N, LBC, RBC, f, g)
shootingIterator::shootingIterator(long double a_in, long double b_in, int N_in, 
                         long double LBC_in, long double RBC_in,
                         long double (*f_in)(long double, long double, long double),
                         long double (*g_in) (long double, long double, long double,
                                              long double*, long double*,
                                              long double*, long double*)) {

    //Most importantly, declare no memory has been allocated yet
    allocated = false;
    //Set up passed values
    a = a_in;
    b = b_in;
    N = N_in;
    LBC = LBC_in;
    RBC = RBC_in;
    f = f_in;
    g = g_in; 
    //Allocate memory
    reallocate();
    //Calculate other needed values
    h = (b - a) / N; 
    t = (RBC - LBC) / (b - a);
    generateXvalues();
    iteration = 1;
    isFirstIteration = true;
    //Set default settings and initial flags
    outputPrecision = 6; 
    useEulers();
    useNewtons();
    streaming = false;
    fileName = "";
    //Push initial conditions for y_system and z_system
    y_approximation[0] = LBC;
    yprime_approximation[0] = t;
    z_approximation[0] = 0.00;
    zprime_approximation[0] = 1.00;
    y_full[0] = LBC;
    yprime_full[0] = t;
    //Record t0
    t_list.push_back(t);
}

//This function clears memory if allocated, then allocates new memory
void shootingIterator::reallocate() {
    if (allocated) deallocate();
    x = new long double [N+1];
    y_approximation = new long double [N+1];
    yprime_approximation = new long double [N+1];
    z_approximation = new long double [N+1];
    zprime_approximation = new long double [N+1];
    y_full = new long double [2*N + 1];
    yprime_full = new long double [2*N + 1];
    allocated = true;
}

//This function clears all heap memory.
void shootingIterator::deallocate() {
    if (allocated) {
        delete[] x;
        delete[] y_approximation;
        delete[] yprime_approximation;
        delete[] z_approximation;
        delete[] zprime_approximation;
        delete[] y_full;
        delete[] yprime_full;
        allocated = false;
    }
}

//This function generates the x-values with the given a, b, and h
void shootingIterator::generateXvalues() {
    int i = -1;
    bool firstOverflow = true;
    do {
        i++;
        if (i < N+1) {
            x[i] = a + i*h;
        }
        else {
            if (firstOverflow) firstOverflow = false; 
                   //ignore first overflow as it often happens (since the
                   //   i++ is at the beginning, which it needs to be for 
                   //   those values of b that are hard to hit precisely)
            else { 
                cout << "Error in generating x values. Not enough room.\n";
                cout << "N = " << N << ", h = " << h << ", \n";
                cout << "Last value added: " << x[N-1] << endl;
                cout << "Value attempted to add: " << (a + i*h) << endl;
            }
        }
    } while (a + i*h <= b);
}

//Opens an output file stream
bool shootingIterator::beginOutput(string outfile) {
    if (streaming) endOutput();
    fileName = outfile;
    fout.open(fileName.c_str());
    if (!(fout)) {
        cout << "Error: file could not be opened for output.\n";
        streaming = false;
    }
    else streaming = true;
    return streaming;
}

//Closes any open file 
bool shootingIterator::endOutput() {
    fout.close();
    streaming = false;
    return true;
}

//Prints the current iteration to the console
void shootingIterator::print() const { 
    string output = getPrintOutput();
    cout << output;
}

//Prints the current iteration to file
void shootingIterator::printToFile() {
    if (streaming) {
        string output = getPrintOutput();
        fout << output;
    }
    else {
        cout << "Error: cannot print to file. No file currently open.\n";
    }
}

//Adds a string to the current file output,
//  adding a blank line before and after
void shootingIterator::addToFile(string out) {
    if (streaming)
        fout << endl << out << endl << endl;
}

//Constructs output for printing the current iteration
string shootingIterator::getPrintOutput() const {
    stringstream output;
    output << endl;
    output << "********** Iteration " << iteration << " **********\n";
    output << fixed << setprecision(outputPrecision);
    output << "t = " << t << "\n";
    output << "\n";
    for (int i = 0; i < N+1; i++) {
        output << "    ";
        output << "x = " << x[i] << " >>\t ";
        output << "y ~=  " << y_approximation[i] << " \t";
        output << "y' ~= " << yprime_approximation[i] << " \t";
        output << "z ~=  " << z_approximation[i] << " \t";
        output << "z' ~= " << zprime_approximation[i] << endl;
    }
    output << "Distance to target = " << error() << endl;
    output << endl;
    return output.str();
}

//Prints a summary of the current iteration to the console
void shootingIterator::printSummary() const {
    stringstream output;
    output << "Iteration " << iteration << ": ";
    output << fixed << setprecision(outputPrecision);
    output << "t = " << t << "; y(b) ~= " << y_approximation[N];
    output << "; target error = " << error() << endl;
    cout << output.str();
}

//Retrieves a single data point from the data
long double shootingIterator::getDataPoint(int column, int i) const {
    long double result = 0.00;
    if (allocated) {
        if ((i >= 0) && (i < N+1)) {
            switch (column) {
                case 0:
                    result = x[i];
                    break;
                case 1: 
                    result = y_approximation[i];
                    break;
                case 2:
                    result = yprime_approximation[i];
                    break; 
                case 3:
                    result = z_approximation[i];
                    break; 
                case 4: 
                    result = zprime_approximation[i];
                    break;
                default:
                    cout << "Error retrieving data point: invalid column.\n";
            }
        }
        else  cout << "Error retrieving data point: index " << i << " out of range.\n";
    }
    else cout << "Error retrieving data point: no memory allocated\n";
    return result;
}

//Calculates the distance between the final y approximation and the RBC
long double shootingIterator::error() const {
    long double result = 0.00;
    if (allocated) result = RBC - y_approximation[N];
    else cout << "Error calculating target error: no memory allocated\n";
    return result;
}

//This function calculates the next iteration
void shootingIterator::startNextIteration() {

    // If this is not the first iteration, calculate a new t,
    // clear the memory, and set up the new iteration
    if (!(isFirstIteration)) {
        t = tValueIterator();
        reallocate();
        generateXvalues();
        iteration++;
        //Push initial conditions for y_system and z_system
        y_approximation[0] = LBC;
        yprime_approximation[0] = t;
        z_approximation[0] = 0.00;
        zprime_approximation[0] = 1.00;
        y_full[0] = LBC;
        yprime_full[0] = t;
        //Record t
        t_list.push_back(t);
    }
    else isFirstIteration = false; //if it was the first iteration, it 
                                   //    will be so no longer

    //Calculate the y_system meshpoints
    // This is the normal method--if set to runge-kutta, we must do extra
    if (selectedMeshPointFunction != RUNGE_KUTTA) {
        for (int i = 0; i < N; i++) {
            meshPointFunction(x[i], y_approximation[i], yprime_approximation[i],
                              y_approximation[i+1], yprime_approximation[i+1], h,
                              &y_system);
        }
    }
    else { //runge-kutta requires half measures for the y-system, so 
           // populate over twice as many mesh points
           //Note for this one only we must pass an h/2 instead of h as 
           // the step size substitution
        for (int i = 0; i < 2*N; i++) {
            meshPointFunction(a + ((h*i)/2), y_full[i], yprime_full[i],
                              y_full[i+1], yprime_full[i+1], (h/2),
                              &y_system);
        }
        //Copy over full measures to the normal array
        for (int i = 0; i < N+1; i++) {
             y_approximation[i] = y_full[i*2];
             yprime_approximation[i] = yprime_full[i*2];
        }
    }

    //Record end y and target error
    error_list.push_back(error());
    endY_list.push_back(getDataPoint(1, N));

    //Calculate the z_system meshpoints
    for (int i = 0; i < N; i++) {
        meshPointFunction(x[i], z_approximation[i], zprime_approximation[i],
                          z_approximation[i+1], zprime_approximation[i+1], h,
                          &z_system);
    }

    //If currently streaming to file, print to file
    if (streaming) printToFile();

}

//Generate a mesh point for a system using runge-kutta method
void shootingIterator::rungeKutta(long double x, long double u1, long double u2,
          long double& u1_next, long double& u2_next, long double hh,
          long double (*system_ptr)(int, long double, long double, long double)) {
    long double k1[2];
    long double k2[2];
    long double k3[2];
    long double k4[2];
    //Note: all the k1's must be computed before any of the k2's can be computed
    for (int i = 0; i < 2; i++) 
        k1[i] = hh * system_ptr(i, x, u1, u2);
    //All the k1's computed, now we can compute the k2's
    for (int i = 0; i < 2; i++)
        k2[i] = hh * system_ptr(i, x + (hh/2), u1 + 0.5*k1[0], u2 + 0.5*k1[1]);
    //All the k2's computed, now we can compute the k3's
    for (int i = 0; i < 2; i++)
        k3[i] = hh * system_ptr(i, x + (hh/2), u1 + 0.5*k2[0], u2 + 0.5*k2[1]);
    //All the k3's computed, now we can compute the k4's
    for (int i = 0; i < 2; i++)
        k4[i] = hh * system_ptr(i, x + hh, u1 + k3[0], u2 + k3[1]);
    //Now lets construct the meshpoint
    u1_next = u1 + ((k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6.00);
    u2_next = u2 + ((k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6.00);
}

//Generate a mesh point for a system using eulers method
void shootingIterator::eulers(long double x, long double u1, long double u2,
          long double& u1_next, long double& u2_next,
          long double (*system_ptr)(int, long double, long double, long double)) {
    u1_next = u1 + (h * system_ptr(0, x, u1, u2));
    u2_next = u2 + (h * system_ptr(1, x, u1, u2));
}

//Use newton's method to find the next t-value
long double shootingIterator::newtonsMethod() {
    long double result = t - ((y_approximation[N] - RBC) / z_approximation[N]);
    return result;
}

/**
 * The following functions evaluate the given user-defined system
 * which begins from the BVP involving  y"(x,t) = F(x, y, y').
 * The system of two initial value problems is constructed by
 * letting y = y1 and y' = y2, z = z1 and z' = z2, where  
 * z"(x,t) = df/dy * z(x,t) + df/dy' * z'(x,t)  (partial deriv's)
 * Thus the system expected is in the form:
 *
 *    y1' = y2
 *    y2' = f(y1, y2, x)
 *
 *    z1' = z2
 *    z2' = g(z1, z2, x)
 *
 * where f and g are functions defined by the user in the main file
 */
//Internal functions that reference f and g only as needed, to simplify for user
long double shootingIterator::y_system(int index, long double x, long double y1, long double y2) {
    long double result = 0.00;
    switch (index) {
        case 0: 
            result = y2;
            break;
        case 1:
            result = f(x, y1, y2);
            break;
    }
    return result;
}
long double shootingIterator::z_system(int index, long double x, long double z1, long double z2) {
    long double result = 0.00;
    switch (index) {
        case 0:
            result = z2;
            break;
        case 1:
            result = g(x, z1, z2, y_approximation, yprime_approximation, y_full, yprime_full);
            break;
    }
    return result;
}

//Calls the selected meshpoint generator function
void shootingIterator::meshPointFunction(long double x, long double u1, long double u2,
          long double& u1_next, long double& u2_next, long double hh,
          long double (*system_ptr)(int, long double, long double, long double)) {

    if (selectedMeshPointFunction == RUNGE_KUTTA)
        rungeKutta(x, u1, u2, u1_next, u2_next, hh, system_ptr);

    else if (selectedMeshPointFunction == EULERS_METHOD) 
        eulers(x, u1, u2, u1_next, u2_next, system_ptr);

    else //default
        eulers(x, u1, u2, u1_next, u2_next, system_ptr);

}

//Calls the selected t-value iterator function
long double shootingIterator::tValueIterator() {
    //right now there's only one option
    return newtonsMethod();
}


























