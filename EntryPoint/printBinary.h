#ifndef PRINTBINARY_H_INCLUDED
#define PRINTBINARY_H_INCLUDED



#endif // PRINTBINARY_H_INCLUDED

void printBinary(unsigned short int val)
{
    for(int i = 15; i >= 0; i--) {
        if( val & ( 1 << i )) {
            cout << 1;
        } else {
            cout << 0;
        }
    }
}
