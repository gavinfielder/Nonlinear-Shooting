/*
 * vec.cpp
 * Russell Gavin Fielder, CSU Chico, Math 461 Numerical Analysis
 * Prof. Sergei Fomin
 * Version 4/26/2017
 *
 * This file holds the implementation for the class vec, an n-dimensional numerical vector
 */

#include <iostream>
#include "vec.h"

using namespace std;

/*
//Operator for multiplying a vector by a scalar - causes memory errors
vec operator*(long double s, vec v) {
    vec newVec(v.dim);
    for (int i = 0; i < v.dim; i++)
        newVec[i] = (v.e[i]) * s;
    return newVec;
}
vec operator*(vec v, long double s) {
    vec newVec(v.dim);
    for (int i = 0; i < v.dim; i++)
        newVec[i] = (v.e[i]) * s;
    return newVec;
}
*/

//********************************************
//****** class vec function definitions ******
//********************************************


//Default constructor. Does not allocate memory. Try not to use
vec::vec() {
    memoryAllocated = false;
    dim = 0;
}
//Constructor for objects of type vec
vec::vec(int n) {
    memoryAllocated = false;
    dim = n;
    e = new long double [dim];
    memoryAllocated = true;
    for (int i = 0; i < n; i++)
        e[i] = 0.00;
}
//Overloaded constructor for initializing with a long double []
vec::vec(long double init[], int n) {
    memoryAllocated = false;
    dim = n;
    e = new long double [dim];
    memoryAllocated = true;
    for (int i = 0; i < n; i++) 
        e[i] = init[i];
}
//Destructor
vec::~vec() {
    if (memoryAllocated) {
        delete[] e;
        memoryAllocated = false;
        e = NULL;
    }
}
//Allocates memory, for vecs created with default constructor
void vec::initialize(int n) {
    dim = n;
    if (memoryAllocated) {
        delete[] e;
        memoryAllocated = false;
    }
    e = new long double [dim];
    memoryAllocated = true;
    for (int i = 0; i < n; i++)
        e[i] = 0.00;
}
/* operators current cause memory errors
//Operator for adding two vectors together
vec vec::operator+(vec other) {
    vec newVec(this->dim);
    if (this->dim != other.dim) //mismatched dimensions, return first vector itself
        return newVec; 
    // Avoid corruption
    if ((!(this->hasData())) || (!(other.hasData()))) {
        cout << "Corruption detected in operator +" << endl;
        return newVec;
    }
    for (int i = 0; i < this->dim; i++)
        newVec.e[i] = this->e[i] + other.e[i];
    return newVec;
}
//Operator for subtracting two vectors
vec vec::operator-(vec other) {
    vec newVec(this->dim);
    if (this->dim != other.dim) //mismatched dimensions, return first vector itself
        return newVec; 
    // Avoid corruption
    if ((!(this->hasData())) || (!(other.hasData()))) {
        cout << "Corruption detected in operator -" << endl;
        return newVec;
    }
    for (int i = 0; i < this->dim; i++)
        newVec.e[i] = this->e[i] - other.e[i];
    return newVec;
}
//Operator for assigning a vector
vec& vec::operator=(vec other) {
    //Avoid corruption
    if (!(other.hasData())) {
        cout << "Corruption detected in operator =" << endl;
        return *this;
    }
    //Clear memory
    delete[] e;
    memoryAllocated = false;
    //Resize and reallocate
    dim = other.dim;
    e = new long double [dim];
    memoryAllocated = true;
    //Copy memory
    for (int i = 0; i < dim; i++)
        e[i] = other.e[i];
    return *this;
} */
//Operator for accessing vector components directly
long double& vec::operator[](int index) {
    if ((index >= dim) || (index < 0)) { 
        cout << "Index " << index << " out of range error in operator [] (size is " << dim << "). Returning e[0]." << endl;
        return e[0];
    }
    if (!(hasData())) {
        cout << "error in operator []: requested vector has no memory allocated. Returning e[0]" << endl;
        return e[0];
    }
    else return (e[index]);
}
//Prints the vector
void vec::print() const {
    printm();
}
//Prints the vector in wolfram language
void vec::printm() const {
    if (memoryAllocated) {
        cout << "{";
        int i = 0;
        for (i; i < dim - 1; i++)
           cout << e[i] << ", ";
        cout << e[i] << "}";
    }
    else cout << "Error in printm(): No memory allocated." << endl;
}

//Returns whether memory is currently allocated
bool vec::hasData() {
    return memoryAllocated;
}

//Deallocates memory
void vec::deallocate() {
    delete[] e;
    memoryAllocated = false;
}





//********************************************
//**** class veclist function definitions ****
//********************************************

//Constructor
veclist::veclist(int numberOfVectors, int dimension) {
    if (dimension <= 0) {
        cout << "Error: cannot create veclist of size 0. Setting size to 1." << endl;
        dimension = 1;
    }
    arr = new vec* [numberOfVectors];
    size = numberOfVectors;
    for (int i = 0; i < size; i++) {
        arr[i] = new vec(dimension);
    }
}
//Destructor
veclist::~veclist() {
    if (size > 0) {
        delete[] arr;
        size = 0;
    }
}
//Get vec in list
vec& veclist::item(int index) {
    if ((index >= 0) && (index < size))
        return (*(arr[index]));
    else { 
        cout << "Error in veclist.item(): index out of range" << endl;
        return *(new vec);
    }
}
//Operator for accessing list items directly
vec& veclist::operator[](int index) {
    if ((index >= size) || (index < 0)) { 
        cout << "Index " << index << " out of range error in veclist::operator [] (size is " << size << "). Returning item 0." << endl;
        return *(arr[0]);
    }
    else return *(arr[index]);
}  
//Deallocate all memory and clear list
void veclist::clear() {
    for (int i = 0; i < size; i++) {
        (*(arr[i])).deallocate();
        delete arr[i];
        arr[i] = NULL;
    }
    delete[] arr;
    size = 0;
}













