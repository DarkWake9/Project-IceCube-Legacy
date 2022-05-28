#include <iostream>
#include <time.h>
using  namespace std;

int main()
{
    int nsims = 100000;
    int rounds [nsims];
    float p = 0.5;
    for (int i = 0; i <= nsims; i++)    {
        
        int r = 0;
        int nloss = 0;
        
        while (nloss != 2)  {
        
            r+=1;
            
            if ((float) rand()/RAND_MAX < p){ nloss = 0;}

            else{ nloss += 1;}
        rounds[i] = r;
        
        }
    }

    double m = 0;
    for (int i = 0; i <= nsims; i++)    {
        
        m  = m + rounds[i];
    
    }
    m = m / (float) nsims;

    cout << m <<endl<<"6";
    return 0;
}
