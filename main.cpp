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

static uintmax_t mysqrt(uintmax_t n)
{
    if(n > 1)
    {
        uintmax_t x = 4 * n, lastx = x - 2;
        while(abs((intmax_t)x - (intmax_t)lastx) > 1)
        {
            lastx = x;
            x = (x + n * 4 / x) >> 1;
        }
        return x / 2;
    }
    else
        return n;
}

void test_sqrt()
{
    uintmax_t limit = 0x10000; // tested up to 0x10000000
    for(uintmax_t n = 0; n < limit; n++)
    {
        if((n + 1) % 0x100 == 0)
            cout << setw(15) << n << " " << (100 * n / limit) << "%\x1b[K\r" << flush;
        uintmax_t result1 = (uintmax_t)floor_sqrt((fixed_width_uint<32>)n), result2 = mysqrt(n);
        if(result1 != result2)
        {
            cout << n << "\x1b[K\n " << result1 << "\n " << result2 << endl;
        }
    }
}

template <size_t N>
static fixed_width_uint<N> mydiv(fixed_width_uint<N> n_in, fixed_width_uint<N> d_in)
{
    // change to use newton's method: http://en.wikipedia.org/wiki/Division_algorithm#Pseudocode
    typedef fixed_width_uint<N> uint;
    typedef fixed_width_uint<N * 2> wuint;
    if(d_in == (uint)0)
        return (uint)(1U / 0U); // division by zero on purpose
    if(d_in == (uint)1)
        return n_in;
    if(d_in == (uint)2)
        return n_in >> 1;
    if(d_in > n_in)
        return (uint)0;
    if((wuint)d_in << 1 > (wuint)n_in)
        return (uint)1;
    intmax_t extra_fraction_bits = 2;
    intmax_t fraction_bits = ilog2(d_in) + 1 + extra_fraction_bits;
    wuint n = (wuint)n_in << extra_fraction_bits, d = (wuint)d_in << extra_fraction_bits;
    wuint target = (wuint)1 << fraction_bits;
    wuint two = (wuint)2 << fraction_bits;
    while(target + (wuint)1 < d || target - (wuint)1 > d)
    {
        wuint factor = two - d;
        n = n * factor >> fraction_bits;
        d = d * factor >> fraction_bits;
    }
    wuint q = n >> fraction_bits;
    q++;
    if((wuint)n_in < q * (wuint)d_in)
        q--;
    if((wuint)n_in < q * (wuint)d_in)
        q--;
    return (uint)q;
}

void test_div()
{
    typedef fixed_width_uint<32> uint;
#define SINGLE_TEST_VALUE 0x641003
    uintmax_t limit = 0x1000000;
#ifndef SINGLE_TEST_VALUE
    for(uintmax_t n = 0; n < limit; n++)
    {
#else
    do
    {
        constexpr uintmax_t n = SINGLE_TEST_VALUE;
#endif
        if((n + 1) % 0x10000 == 0)
            cout << setw(15) << hex << n << dec << " " << (100 * n / limit) << "%\x1b[K\r" << flush;
        uint num = (uint)(n >> 12);
        uint d = (uint)(n & 0xFFF);
        if(d == (uint)0)
            continue;
        uint result1 = num / d, result2 = mydiv(num, d);
        if(result1 != result2)
        {
            cout << hex << n << dec << "\x1b[K\n " << result1 << "\n " << result2 << endl;
        }
#ifndef SINGLE_TEST_VALUE
    }
#else
#undef SINGLE_TEST_VALUE
    }
    while(false);
#endif
}

int main(int argc, char **argv)
{
    typedef ieee754_soft_float_64 fp_type;

    auto pow10v = (fp_type)(uintmax_t)1;

    for(int i = 1; i < std::numeric_limits<fp_type>::digits10; i++)
        pow10v *= (fp_type)(uintmax_t)10;

    dumpf(fp_type::pi());
    dumpbi((fp_type::word_type)(fp_type::pi() * pow10v));
    auto v = fp_type::log10_2();
    dumpf(v);
    dumpbi((fp_type::word_type)(v * pow10v));

    test_div();

    return 0;
}
