#ifndef IR_INSTRUCTION_H_INCLUDED
#define IR_INSTRUCTION_H_INCLUDED

#include <string>
#include <vector>

template <typename T>
struct register_or_literal
{
    typedef T RegisterType;
    bool isLiteral;

};

template <typename T>
struct ir_instruction
{
    typedef T RegisterType;
    virtual string toString() const;
};

#endif // IR_INSTRUCTION_H_INCLUDED
