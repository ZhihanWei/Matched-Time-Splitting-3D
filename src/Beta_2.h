#ifndef __Beta_2_H__
#define __Beta_2_H__

#include "Constant.h"
#include "Beta.h"

class Beta_2: public Beta
{
public:
    Beta_2();
    
    virtual double Inside(Doub_I, Doub_I, Doub_I) const;
    virtual double Outside(Doub_I, Doub_I, Doub_I) const;
};

#endif

