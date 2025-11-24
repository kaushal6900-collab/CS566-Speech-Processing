#include "stdafx.h"
#include "stdio.h"
#include "process.h"


const long long tnrg = 20000; //i have assumed(trial/error) it by observing the wave amplitude

int _tmain(int argc, _TCHAR* argv[])
{
    int x, i = 0;
    int arr[60001] = {0}; 
    int a[60001] = {0};
    int b[60001] = {0};
    int size;
    int err,c=0;

    FILE *f1 = NULL;
    
    err = fopen_s(&f1, "ysn.txt", "r");
    
    if (err != 0) {
        printf("\n Cannot Open file 'ysn.txt'");
        exit(1);
    }

    int p=0,q=0;
    while (fscanf_s(f1, "%d", &x) == 1) {
        if(i < 60000){
            arr[i] = x;
            
            long long energy = (long long)x * x;

            
            if (c == 0) {
                if (energy > tnrg) {///end of initial noise and start of 1st word
                    
                    a[p] = x;
                    p++;
                } else {
                    
                    if (p > 0) {
                        c = 1; //end of 1st word
                    }
                    
                }
            } else if (c == 1) {
                if (energy > tnrg) {//start of second word
                    
                    b[q] = x;
                    q++;
                } else {
                    
                    if (q > 0) {
                        break; //sice i have recorded only two words so after end of second word will stop storing data
                    }
                   
                }
            }
            i++;
        } else {
            break; 
        }
    }

    int end1 = p;
    int end2 = q;
    fclose(f1);
    size = i;
    
    int zcr1 = 0, zcr2 = 0;
    
    // ZCR calculation for 1st word
    for(int k = 0; k < end1 - 1; k++){
        if ((a[k] >= 0 && a[k+1] < 0) || (a[k+1] >= 0 && a[k] < 0)) {
            zcr1++;
        }
    }

    // ZCR calculation for 2nd word
    for(int k = 0; k < end2 - 1; k++){
        if ((b[k] >= 0 && b[k+1] < 0) || (b[k+1] >= 0 && b[k] < 0)) {
            zcr2++;
        }
    }
    
    
    
    if (zcr1 > zcr2) {
        printf("Result: First word is 'yes'\n");
    }
    else {
        printf("Result: Second word is 'yes'\n");
    }
    
    return 0;
}