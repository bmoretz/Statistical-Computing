#include <iostream>
#include <string>

#include "storageSize.h"
#include "printBinary.h"

using namespace std;

int main()
{
    string value;

    cout << "Enter a number to convert to binary: ";
    cin >> value;

    printBinary(stoi(value));

    string in;
    cin >> in;

    return 0;
}
