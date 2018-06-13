/*
 * vec.h
 * Russell Gavin Fielder, CSU Chico, Math 461 Numerical Analysis
 * Prof. Sergei Fomin
 * Version 4/24/2017
 *
 * This file is the header for class vec, an n-dimensional numerical vector
 */

#ifndef VEC_H
#define VEC_H

using namespace std;


//Objects of this class will store n-dimensional numerical vectors
class vec {
    public:
        // Constructors
        vec();                              // default constructor, does not allocated memory. Try not to use.
        vec(int n);                         // sets up a 0 vector of n dimensions
        vec(long double [], int n);         // initializes the vector with an array of long doubles
        // Destructor
        ~vec();
        // Member variables
        int dim;                            // stores the number of dimensions
        long double * e;                    // stores the components of the vector
        // Member functions
        void print() const;                 // prints the vector; right now this just calls printm()
        void printm() const;                // prints the vector in wolfram language
        bool hasData();                     // Accessor for memoryAllocated
        void initialize(int n);             // Allocates memory (for vectors created with default constructor)
        void deallocate();                  // Deallocates memory
        // Operators
        //vec operator+(vec);                 // adds two vectors - causes memory errors
        //vec& operator=(vec);                // assigns a vector - causes memory errors
        //vec operator-(vec);                 // subtracts two vectors - causes memory errors
        long double& operator[](int index); // retrieves a vector component directly
    private:
        bool memoryAllocated; 
};
 
class veclist {
    public:
        veclist(int numberOfVectors, int dimension);
        ~veclist();
        vec& item(int);
        vec& operator[](int index);         // retrieves a vector from the list
        void clear();                       // deallocates all memory
        int getSize() { return size; }
    private:
        vec** arr;
        int size;
};


//Operators for scalar multiplication
//vec operator*(long double, vec); //causes memory errors
//vec operator*(vec, long double); //causes memory errors


#endif












