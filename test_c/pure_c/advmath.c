#include <stdio.h>

int add(int a, int b){
   int c = a + b;
   printf("Hello from C, received a=%d and b=%d and found c=%d.\nLet's wrap this up so we can do something interesting!\n", a, b, c);
   return c;
}

/* int main() {
   // printf() displays the string inside quotation
   printf("Hello, World!");
   return 0;
} */
