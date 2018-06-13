/*
 * shootingIterator.h
 * Russell Gavin Fielder, CSU Chico, Math 461 Numerical Analysis
 * Prof. Sergei Fomin
 * Version 4/30/2017
 *
 * This file is the header for class shootingIterator, which handles the
 * nonlinear shooting method for second-degree boundary value problems
 */

#ifndef SHOOTER_H
#define SHOOTER_H

#define EULERS_METHOD 0
#define RUNGE_KUTTA 1
#define NEWTONS_METHOD 0
 
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

class shootingIterator {
    private:
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
        static long double y_system(int index, long double x, long double y1, long double y2);
        static long double z_system(int index, long double x, long double z1, long double z2);
        //Note pointers to functions f and g are declared in shootingIterator.cpp and used
        //   by y_system and z_system. 

        //Parallel arrays to store mesh points
        static long double * x;                    //Array to store x values (column 0)
        static long double * y_approximation;      //Array to store y approximations (col 1)
        static long double * yprime_approximation; //Array to store y' approximations (col 2)
        static long double * z_approximation;      //Array to store z approximations (col 3)
        static long double * zprime_approximation; //Array to store z' approximations (col 4)
        static long double * y_full;
        static long double * yprime_full;

        //BVP constants
        long double a, b, LBC, RBC, h;
        int N;

        /**
         * Methods of mesh point approximation; these functions take a mesh point
         * (u1, u2) and generate a new mesh point of the passed system. 
         * The new mesh point is stored in (u1_next, u2_next).
         *     simple form:  rungeKutta(x, u1, u2, &u1_next, &u2_next, y_system)
        */
        void rungeKutta(long double x, long double u1, long double u2,
                 long double& u1_next, long double& u2_next, long double hh,
                 long double (*system_ptr)(int, long double, long double, long double));
        void eulers(long double x, long double u1, long double u2,
                 long double& u1_next, long double& u2_next,
                 long double (*system_ptr)(int, long double, long double, long double));

        //Methods of t-value approximation; these functions evaluate t and the
        //   current iteration's data to generate a new t-value
        long double newtonsMethod();

        //Wrapper functions and program settings for method of approximation
        void meshPointFunction(long double, long double, long double,
                 long double&, long double&, long double,
                 long double (*system_ptr)(int, long double, long double, long double));
        long double tValueIterator();
        static int selectedMeshPointFunction;
        static int selectedTValueIterator;

        //Iteration data
        long double t;                   // Shooting method parameter t
        vector<long double> t_list;      // Tracks history of t values
        vector<long double> endY_list;   // Tracks history of end y values
        vector<long double> error_list;  // Tracks history of (RBC - endY) target errors
        int iteration;                   // Stores the current iteration number

        //Output file handling
        ofstream fout;
        string fileName;    
        bool streaming;     

        //Standard output handling (both file and console)
        string getPrintOutput() const; // Generates standard output
        int outputPrecision;           // Number of decimal places to display

        //Other
        bool allocated;         // Internal flag for whether memory is currently allocated
        void generateXvalues(); // Generates a + i*h values, populating the x[] array
        bool isFirstIteration;  // Flag for first iteration run

    public:
        //Constructor
        // simple form: shootingIterator(a, b, N, LBC, RBC, f, g)
        shootingIterator(long double a_in, long double b_in, int N_in, 
                         long double LBC_in, long double RBC_in,
                         long double (*f_in)(long double, long double, long double),
                         long double (*g_in) (long double, long double, long double,
                                              long double*, long double*,
                                              long double*, long double*));

        //Functions to switch methods of approximation
        void useEulers() { selectedMeshPointFunction = EULERS_METHOD; }
        void useRungeKutta() { selectedMeshPointFunction = RUNGE_KUTTA; }
        void useNewtons() { selectedTValueIterator = NEWTONS_METHOD; }

        //Functions to control output stream
        bool beginOutput(string outfile);   // opens a file for output
        bool endOutput();                   // closes any open file 
        void print() const;                 // prints the current iteration
        void printToFile();                 // prints the iteration to file; note it
                                            //    does automatically when streaming
        void addToFile(string out);         // prints an additional string to file
        void setPrecision(int prec) { outputPrecision = prec; }
        void printSummary() const;          // prints a summary of the iteration

        //Information access
        int iterationNumber() const { return iteration; }
        string getFileName() const { return fileName; }
        int numMeshPoints() const { return N; }
        long double stepSize() const { return h; }
        long double getDataPoint(int column, int i) const;
        long double error() const;
        long double getT() const { return t; }
 
        //Control functions
        void startNextIteration();

        //Memory management
        void reallocate();
        void deallocate();
        bool isAllocated() const { return allocated; }

};



#endif









