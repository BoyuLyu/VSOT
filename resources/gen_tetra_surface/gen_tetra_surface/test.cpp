#define _USE_MATH_DEFINES
#include "stdio.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
using namespace std;
void main() {
    int* xxshift = new int[192];
    int* yyshift = new int[192];
    int* zzshift = new int[192];
    for (size_t i = 0; i < 192; i++)
    {
        xxshift[i] = 0;
        yyshift[i] = 0;
        zzshift[i] = 0;

    }
    int count = 0;
    for (int i = -1; i < 2; i++)
    {
        for (int j = -1; j < 2; j++)
        {
            for (int k = -1; k < 2; k++)
            {
                if (i == 0 && j == 0 && k != 0)
                {
                    for (int m = -1; m < 1; m++)
                    {
                        for (int n = -1; n < 1; n++)
                        {
                            xxshift[count * 4] = 0;
                            xxshift[count * 4 + 1] = 0;
                            xxshift[count * 4 + 2] = m;
                            xxshift[count * 4 + 3] = 0;

                            yyshift[count * 4] = 0;
                            yyshift[count * 4 + 1] = 0;
                            yyshift[count * 4 + 2] = n;
                            yyshift[count * 4 + 3] = n;

                            zzshift[count * 4] = 0;
                            zzshift[count * 4 + 1] = 2*k;
                            zzshift[count * 4 + 2] = k;
                            zzshift[count * 4 + 3] = k;

                            count++;
                            xxshift[count * 4] = 0;
                            xxshift[count * 4 + 1] = 0;
                            xxshift[count * 4 + 2] = m;
                            xxshift[count * 4 + 3] = m;

                            yyshift[count * 4] = 0;
                            yyshift[count * 4 + 1] = 0;
                            yyshift[count * 4 + 2] = n;
                            yyshift[count * 4 + 3] = 0;

                            zzshift[count * 4] = 0;
                            zzshift[count * 4 + 1] = 2 * k;
                            zzshift[count * 4 + 2] = k;
                            zzshift[count * 4 + 3] = k;

                            count++;
                        }
                    }

                }
                else if (i==0 && j !=0 && k ==0)
                {
                    for (int m = -1; m < 1; m++)
                    {
                        for (int n = -1; n < 1; n++)
                        {
                            xxshift[count * 4] = 0;
                            xxshift[count * 4 + 1] = 0;
                            xxshift[count * 4 + 2] = m;
                            xxshift[count * 4 + 3] = 0;

                            yyshift[count * 4] = 0;
                            yyshift[count * 4 + 1] = 2*j;
                            yyshift[count * 4 + 2] = j;
                            yyshift[count * 4 + 3] = j;

                            zzshift[count * 4] = 0;
                            zzshift[count * 4 + 1] = 0;
                            zzshift[count * 4 + 2] = n;
                            zzshift[count * 4 + 3] = n;

                            count++;
                            xxshift[count * 4] = 0;
                            xxshift[count * 4 + 1] = 0;
                            xxshift[count * 4 + 2] = m;
                            xxshift[count * 4 + 3] = m;

                            yyshift[count * 4] = 0;
                            yyshift[count * 4 + 1] = 2*j;
                            yyshift[count * 4 + 2] = j;
                            yyshift[count * 4 + 3] = j;

                            zzshift[count * 4] = 0;
                            zzshift[count * 4 + 1] = 0;
                            zzshift[count * 4 + 2] = n;
                            zzshift[count * 4 + 3] = 0;
                            count++;
                        }
                    }
                }
                else if (i != 0 && j == 0 && k == 0)
                {
                    for (int m = -1; m < 1; m++)
                    {
                        for (int n = -1; n < 1; n++)
                        {
                            xxshift[count * 4] = 0;
                            xxshift[count * 4 + 1] = 2*i;
                            xxshift[count * 4 + 2] = i;
                            xxshift[count * 4 + 3] = i;

                            yyshift[count * 4] = 0;
                            yyshift[count * 4 + 1] = 0;
                            yyshift[count * 4 + 2] = m;
                            yyshift[count * 4 + 3] = 0;

                            zzshift[count * 4] = 0;
                            zzshift[count * 4 + 1] = 0;
                            zzshift[count * 4 + 2] = n;
                            zzshift[count * 4 + 3] = n;

                            count++;
                            xxshift[count * 4] = 0;
                            xxshift[count * 4 + 1] = 2*i;
                            xxshift[count * 4 + 2] = i;
                            xxshift[count * 4 + 3] = i;

                            yyshift[count * 4] = 0;
                            yyshift[count * 4 + 1] = 0;
                            yyshift[count * 4 + 2] = m;
                            yyshift[count * 4 + 3] = m;

                            zzshift[count * 4] = 0;
                            zzshift[count * 4 + 1] = 0;
                            zzshift[count * 4 + 2] = n;
                            zzshift[count * 4 + 3] = 0;
                            count++;
                        }
                    }
                }

            }
        }
    }
    ofstream myfile;
    myfile.open("xxshift.txt");
    for (size_t i = 0; i < 192; i++)
    {
        myfile << xxshift[i] <<endl;
    }
    myfile.close();
}