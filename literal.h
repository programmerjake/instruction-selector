#ifndef LITERAL_H_INCLUDED
#define LITERAL_H_INCLUDED

#include <cstdint>

struct literal
{
    union ValueType
    {
        double f64;
        float f32;

    };
    ValueType value;
};

#endif // LITERAL_H_INCLUDED
