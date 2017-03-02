#ifndef __Beta_1_H__
#define __Beta_1_H__

#include "Constant.h"
#include "Beta.h"

class Beta_1: public Beta
{
public:
    virtual double Inside(Doub_I, Doub_I, Doub_I) const;
    virtual double Outside(Doub_I, Doub_I, Doub_I) const;
};

#endif

