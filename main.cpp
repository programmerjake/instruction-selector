#include <iostream>
#include <iomanip>
#include "float_types.h"

using namespace std;

template <typename T, T V>
struct template_value final
{
    static constexpr T value = V;
};

int main(int argc, char **argv)
{
    cout << hex;
    auto p = (ieee754_soft_float_32((long double)1.)).get_bits();
    cout << "0x" << setw(decltype(p)::bit_count / 4) << setfill('0') << p << endl;
    return 0;
}
