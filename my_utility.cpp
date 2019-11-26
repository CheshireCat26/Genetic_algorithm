//
// Created by cheshirecat on 11/15/19.
//

#include "my_utility.h"

unsigned int factorial(unsigned int n)
{
    unsigned int ret = n;
    n--;
    for (; n >= 1; n--)
        ret *= n;
    return ret;
}
