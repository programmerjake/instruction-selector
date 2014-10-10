#ifndef INSTRUCTION_LIST_H_INCLUDED
#define INSTRUCTION_LIST_H_INCLUDED

#include <memory>
#include <vector>

template <typename T>
struct basic_blocks;

template <typename T>
struct basic_block
{
    std::vector<T> instructions;
    std::vector<std::weak_ptr<basic_block<T>>> sourceBlocks;
    std::vector<std::weak_ptr<basic_block<T>>> targetBlocks;
    std::weak_ptr<basic_blocks<T>> allBlocks;
};

template <typename T>
struct basic_blocks
{
    std::vector<std::shared_ptr<basic_block<T>>> blocks;
};

#endif // INSTRUCTION_LIST_H_INCLUDED