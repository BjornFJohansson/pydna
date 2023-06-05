#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int smallest_rotation_index(char *s) {
    int n = strlen(s);
    int i = 0, j = 1, k = 0;
    while (i < n && j < n && k < n) {
        int cmp = s[(i+k)%n] - s[(j+k)%n];
        if (cmp == 0) {
            k++;
        } else {
            if (cmp > 0) {
                i = i + k + 1;
            } else {
                j = j + k + 1;
            }
            if (i == j) {
                j++;
            }
            k = 0;
        }
    }
    return i < j ? i : j;
}

void reverse(char *a, int n) {

  for (int i=0, j=n-1; i < j; i++, j--) {
      a[i] = a[i] ^ a[j];
      a[j] = a[j] ^ a[i];
      a[i] = a[i] ^ a[j];
  }
}


char * smallest_rotation(char *myString)
{
    int n = strlen(myString);
    int k = n - smallest_rotation_index(myString);

    reverse(&myString[n-k], k); // reverse end of array
    reverse(myString, n-k);     // reverse beginning of array
    reverse(myString, n);       // reverse entire array

    return myString;

}

int main() {
    char s[] = "cabcdab";

    char * s2 = smallest_rotation(s);

    printf("%s\n", s2);

    return 0;
}
