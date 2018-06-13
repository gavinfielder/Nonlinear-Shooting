/*
 * 12_NonlinearShooting.cpp
 * Russell Gavin Fielder, CSU Chico, Math 461 Numerical Analysis
 * Prof. Sergei Fomin
 * Version 4/30/2017
 *
 * This program applies the Shooting Method to a second order nonlinear BVP
 */

/** 
 *         Instructions for Use
 *
 *  1. Ensure files "shootingIterator.h", "shootingIterator.cpp" present
 *  2. Enter BVP Constants and program settings starting at line 35
 *  3. Define the 2nd order differential equation starting at line 158
 *  4. Compile. (For g++, the only command is "make")
 *  5. Run the output ("./nonlinearshooting"). 
*/

//****************** do not change ******************

#include <iostream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include "shootingIterator.h"

using namespace std;

//Misc constants
const long double PI = 3.14159265358979323846; 


//********************************************************************
//************************ Program Settings **************************
//********************************************************************

//X-Interval [a,b]
const long double a = 0.00;
const long double b = PI;
//Left and right dirichlet boundary conditions
const long double LBC = 2.00;      
const long double RBC = 2.00;     

//Number of subintervals
const int N = 20;
long double h; // h will be calculated automatically

//Method of approximation - if false, use eulers
const bool USE_RUNGE_KUTTA = true; 

//Standard output precision used in console output as well as
const int PRECISION = 15;
const int MAX_ITERATIONS = 100;

//Desired shot accuracy 
const long double TOLERANCE = 0.0001; 

//Set this to true for meshpoint error analysis option 
//   Provide exact solution below
const bool SOLUTION_KNOWN = true;  
const bool SUPPRESS_CONSOLE_MESHPOINT_ERROR_OUTPUT = true;

//Meshpoint error analysis settings
const int ANALYSIS_X_PREC = 7;     //output precision to use for x
const int ANALYSIS_Y_PREC = 12;    //output precision to use for y, w

//********************************************************************
//********************* end program settings *************************
//********************************************************************



//****************** do not change ******************

//Function prototypes
long double f(long double x, long double y, long double yprime);
long double g(long double x, long double y, long double yprime,
              long double * y_approx, long double * yprime_approx,
              long double * y_full, long double * yprime_full);
long double analyzeResults(shootingIterator&, int precX, int precY, 
               bool suppress = SUPPRESS_CONSOLE_MESHPOINT_ERROR_OUTPUT);
long double solution(long double x);
long double dfdy(long double x, long double y, long double yprime);
long double dfdyprime(long double x, long double y, long double yprime);


//Main function
int main() {
    char userChar = 'n';

    //Declare the shooting method handler (see shootingIterator.h)
    shootingIterator shooter(a, b, N, LBC, RBC, &f, &g);
    h = shooter.stepSize();

    //Set program settings
    if (USE_RUNGE_KUTTA) shooter.useRungeKutta();
    else shooter.useEulers();
    shooter.useNewtons();
    shooter.setPrecision(PRECISION);

    //Begin data output to file
    shooter.beginOutput("output.txt");

    //Welcome message
    cout << endl;
    stringstream output;
    output << "*************************************************\n";
    output << "****** Nonlinear Shooting Method results: *******\n";
    output << "*************************************************\n\n";
    output << "Approximating the solution of second-order BVP\n";
    output << "with y(" << a << ") = " << LBC << ", y(" << b << ") = " << RBC << ", with N = " << N << "\n";
    output << "subintervals,  " << a << " <= x <= " << b << " :\n";
    if (USE_RUNGE_KUTTA) output << "Using Runge-Kutta method of order 4 and ";
    else output << "Using Euler's method and\n";
    output << "\nNewton's Method for iterative parameter t\n";
    output << "Target tolerance is " << TOLERANCE << "\n";
    output << "Calculations are being streamed to: " << shooter.getFileName() << ".\n\n";
    cout << output.str();
    shooter.addToFile(output.str());
    cout << "Basic information:\n\n";

    //Start iterating
    do {
        shooter.startNextIteration();
        shooter.printSummary();
    } while ((abs(shooter.error()) > TOLERANCE) 
            && (shooter.iterationNumber() < MAX_ITERATIONS));

    //Done iterating, figure out why
    if (shooter.iterationNumber() >= MAX_ITERATIONS) {
        cout << "\nMaximum number of iterations reached. Recommend switching\n";
        cout << "  approximation methods for faster convergence.\n";
    }
    else {
        //Tolerance reached
        shooter.addToFile("Tolerance achieved. Stopping Iterations.");
        cout << "\nReached desired tolerance at " << shooter.iterationNumber();
        cout << " iterations. ";
        //Check if the user provided a solution to check against
        if (SOLUTION_KNOWN) {
            cout << "Exact solution detected; analyzing meshpoint errors.\n";
            long double maxError = analyzeResults(shooter,
                                     ANALYSIS_X_PREC, ANALYSIS_Y_PREC);
        }
        else cout << endl;
    }
    cout << endl << endl;

    shooter.endOutput();
    shooter.deallocate();

    return 0;
}


//********************************************************************
//************ Define 2nd Order Differential Equation ****************
//********************************************************************

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
 * To use this system, simply input f, df/dy, and df/dy' below:
 * If the actual solution is known, you may enter it into solution(x)
*/

//Here is f(y1, y2, x)
long double f(long double x, long double y, long double yprime) {

    //Edit starting here v
    long double result = (1 - pow(yprime, 2) - y * sin(x)) / 2;
    //End editing        ^
    return result;
}
//Here is the partial derivative of f with respect to y
long double dfdy(long double x, long double y, long double yprime) {

    //Edit starting here v
    long double result = sin(x) / -2.00;
    //End editing        ^
    return result;
}
//Here is the partial derivative of f with respect to y'
long double dfdyprime(long double x, long double y, long double yprime) {
 
    //Edit starting here v
    long double result = yprime * -1;
    //End editing
    return result;
}

//This is the actual solution, if known. This allows checking meshpoint error
long double solution(long double x) {

    //Edit starting here v
    long double result = 2 + sin(x);
    //End editing        ^
    return result;
}


//********************************************************************
//********** end define 2nd order differential equation **************
//********************************************************************


//************************* do not change ****************************

//Here is g(x, z1, z2)
long double g(long double x, long double z, long double zprime,
              long double * y_approx, long double * yprime_approx,
              long double * y_full, long double * yprime_full) {

    //Get approximation mesh points from y-system
    int i = 0;
    long double y, yprime;
    if (!(USE_RUNGE_KUTTA)) {
        i = lround((x-a)/h);
        y = y_approx[i];
        yprime = yprime_approx[i];
    }
    else { //Runge-kutta used. Need to determine whether requested
           //  x-value is a half-measure in approximation arrays
        if (lround(((2*(x-a))/h)) % 2 == 0) {
            //not a half measure, use y_approx
            i = lround((x-a)/h);
            y = y_approx[i];
            yprime = yprime_approx[i];
        }
        else {
            //half measure, get data from full array
            i = lround((((x-a)/h) * 2));
            y = y_full[i];
            yprime = yprime_full[i];
        }
    }

    long double result = dfdy(x, y, yprime) * z + dfdyprime(x, y, yprime) * zprime;
    return result;
}


//Calculates and prints errors for each meshpoint 
long double analyzeResults(shootingIterator& shooted, int x_precision,
                    int y_precision, bool suppressConsoleOutput) {
    stringstream output;
    long double error = 0.00;
    long double maxError = 0.00;
    output << "***************************************************\n";
    output << "*************** Mesh point errors *****************\n";
    output << "***************************************************\n";
    output << endl;
    output << "Iteration number: " << shooted.iterationNumber() << endl;
    output << "t = " << fixed << setprecision(y_precision) << shooted.getT() << endl << endl;
    for (int i = 0; i < N + 1; i++) {
        output << "  x = " << fixed << setprecision(x_precision) << shooted.getDataPoint(0, i);
        output << ", \tw = " << fixed << setprecision(y_precision) << shooted.getDataPoint(1, i);
        output << ", \ty = " << fixed << setprecision(y_precision) << solution(shooted.getDataPoint(0,i));
        error = shooted.getDataPoint(1, i) - solution(shooted.getDataPoint(0,i));
        output << ", \tError = " << scientific << error << endl;
        if (abs(error) > maxError) maxError = abs(error);
    }
    output << endl;

    output << "Maximum absolute meshpoint error is " << scientific << maxError << "\n";  

    //Add this to the file if streaming
    shooted.addToFile(output.str());

    //Print analysis to console
    if (!(suppressConsoleOutput)) {
        cout << endl << output.str();
    }
    else {
        cout << "Maximum absolute meshpoint error is " << scientific << maxError << "\n";
        cout << "Meshpoint error analysis provided with " << shooted.getFileName();
        cout << endl;
    }
  
    return maxError;
}













