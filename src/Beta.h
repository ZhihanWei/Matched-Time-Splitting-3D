#ifndef __Beta_H__
#define __Beta_H__

#include "Constant.h"

class Beta
{
public:
    virtual double Inside(Doub_I, Doub_I, Doub_I) const = 0;
    virtual double Outside(Doub_I, Doub_I, Doub_I) const = 0;
};

#endif
