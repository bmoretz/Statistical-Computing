#ifndef STORAGESIZE_H_INCLUDED
#define STORAGESIZE_H_INCLUDED

#include <iostream>

#endif // STORAGESIZE_H_INCLUDED

using namespace std;

void display_storage()
{
    cout << "The storage allocated for a char is " << sizeof(char)
        << " bytes" << endl;
    cout << "The storage allocated for an unsigned short integer is " << sizeof(unsigned short int)
        << " bytes" << endl;
    cout << "The storage allocated for an integer is " << sizeof(int)
        << " bytes" << endl;
    cout << "The storage allocated for a long integer is " << sizeof(long int)
        << " bytes" << endl;
    cout << "The storage allocated for an unsigned long integer is " << sizeof(unsigned long int)
        << " bytes" << endl;
    cout << "The storage allocated for a float is " << sizeof(float)
        << " bytes" << endl;
    cout << "The storage allocated for a double is " << sizeof(double)
        << " bytes" << endl;
    cout << "The storage allocated for a long double is " << sizeof(long double)
        << " bytes" << endl;
}
