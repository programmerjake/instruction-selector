#include <iostream>
#include <iomanip>
#include "float_types.h"

using namespace std;

template <typename T, T V>
struct template_value final
{
    static constexpr T value = V;
};

template <size_t N>
struct to_float;

template <>
struct to_float<32>
{
    union
    {
        float value;
        uint32_t ivalue;
    };
    to_float(uint32_t v)
        : ivalue(v)
    {
    }
};

template <>
struct to_float<64>
{
    union
    {
        double value;
        uint64_t ivalue;
    };
    to_float(uint64_t v)
        : ivalue(v)
    {
    }
};

template <size_t N>
struct to_word;

template <>
struct to_word<32>
{
    union
    {
        float fvalue;
        uint32_t value;
    };
    to_word(float v)
        : fvalue(v)
    {
    }
};

template <>
struct to_word<64>
{
    union
    {
        double fvalue;
        uint64_t value;
    };
    to_word(double v)
        : fvalue(v)
    {
    }
};

template <size_t N>
void dump_fn(const char *value_name, uintmax_t v)
{
    cout << value_name << ": ";
    cout << "0x" << uppercase << setw(N / 4) << setfill('0') << hex << v << " " << dec;
    //cout << (intmax_t)v << " ";
    cout << setprecision(6) << showpoint << to_float<64>((uintmax_t)((ieee754_soft_float_64)ieee754_soft_float_std<N>::type::from_bits((typename ieee754_soft_float_std<N>::type::word_type)v)).get_bits()).value << " ";
    //cout << v * std::pow(0.5, ieee754_soft_float_std<N>::type::radix_point_shift) << " ";
    cout << endl;
}

template <size_t N>
void dump_fn(const char *value_name, typename ieee754_soft_float_std<N>::type v)
{
    cout << value_name << ": ";
    cout << "0x" << uppercase << setw(N / 4) << setfill('0') << hex << v.get_bits() << " " << dec;
    cout << setprecision(15) << showpoint << to_float<64>((uintmax_t)((ieee754_soft_float_64)v).get_bits()).value << " ";
    cout << endl;
}

template <size_t N>
void dump_fn(const char *value_name, typename ieee754_soft_float_std<N>::type::word_type v)
{
    cout << value_name << ": ";
    cout << "0x" << uppercase << setw(N / 4) << setfill('0') << hex << v << " " << dec;
    cout << v << " ";
    cout << endl;
}

template <size_t N>
void dump_fn(const char *value_name, typename ieee754_soft_float_std<N>::type::signed_type v)
{
    dump_fn<N>(value_name, (typename ieee754_soft_float_std<N>::type::word_type)v);
}

#define dumpf(v) dump_fn<decltype(v)::word_type::bit_count>(#v, v)
#define dumpbi(v) dump_fn<decltype(v)::bit_count>(#v, v)

int main(int argc, char **argv)
{
    typedef ieee754_soft_float_512 fp_type;

    auto pow10v = (fp_type)(uintmax_t)1;

    for(int i = 1; i < std::numeric_limits<fp_type>::digits10; i++)
        pow10v *= (fp_type)(uintmax_t)10;

    dumpf(fp_type::pi());
    dumpbi((fp_type::word_type)(fp_type::pi() * pow10v));

    return 0;
}
