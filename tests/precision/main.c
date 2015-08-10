#include <stdio.h>

int main() {
    double x, y; 
    x = 1.1234567890123456;
    y = 1234567890.123456789012345678901234567890;

    printf("%50e %50f\n", x, x);
    printf("%50e %50f\n", y, y);
    printf("%.50e %.50f\n", x, x);
    printf("%.50e %.50f\n", y, y);

    printf("%.16e %.16f\n", x, x);
    printf("%.16e %.16f\n", y, y);

    printf("%.17e %.17f\n", x, x);
    printf("%.17e %.17f\n", y, y);
    
    return(0);
}
