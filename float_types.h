#ifndef FLOAT_TYPES_H_INCLUDED
#define FLOAT_TYPES_H_INCLUDED

#include <cmath>
#include "int_types.h"

#warning fix ieee754_soft_float

template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize = 0>
class ieee754_soft_float
{
public:
    static constexpr size_t mantissa_bit_count = mantissaBitCount;
    static constexpr size_t exponent_bit_count = exponentBitCount;
    static constexpr bool mantissa_msb_is_implicit = mantissaMsbIsImplicit;
    static constexpr size_t padding_size = paddingSize;
    static constexpr size_t bit_count = mantissa_bit_count + exponent_bit_count + 1 - (mantissa_msb_is_implicit ? 1 : 0);
    typedef fixed_width_uint<bit_count + padding_size> word_type;
    typedef fixed_width_uint<word_type::bit_count * 2> double_word_type;
    typedef fixed_width_int<bit_count + padding_size> signed_type;
    static constexpr size_t sign_bit_shift = bit_count - 1;
    static constexpr size_t exponent_shift = sign_bit_shift - exponent_bit_count;
    static constexpr size_t mantissa_shift = 0;
    static constexpr size_t radix_point_shift = mantissa_bit_count - 1;
    static constexpr word_type sign_bit_mask = (word_type)1 << sign_bit_shift;
    static constexpr word_type exponent_mask = (((word_type)1 << exponent_bit_count) - (word_type)1) << exponent_shift;
    static constexpr word_type exponent_bias = ((word_type)1 << (exponent_bit_count - 1)) - (word_type)1;
    static constexpr word_type inf_nan_biased_exponent = exponent_mask;
    static constexpr word_type denormal_zero_biased_exponent = (word_type)0;
    static constexpr signed_type inf_nan_unbiased_exponent = (signed_type)inf_nan_biased_exponent - (signed_type)exponent_bias;
    static constexpr signed_type denormal_zero_unbiased_exponent = (signed_type)denormal_zero_biased_exponent - (signed_type)exponent_bias;
    static constexpr signed_type denormal_actual_exponent = denormal_zero_unbiased_exponent + (signed_type)1;
    static constexpr signed_type min_denormal_normalized_exponent = denormal_actual_exponent - (signed_type)(intmax_t)mantissa_bit_count;
    static constexpr word_type mantissa_mask = ((word_type)1 << (mantissa_bit_count - (mantissa_msb_is_implicit ? 1 : 0))) - (word_type)1;
    static constexpr word_type mantissa_implicit_bit = mantissa_msb_is_implicit ? (word_type)1 << (mantissa_bit_count - 1) : (word_type)0;
    static constexpr word_type nan_word = exponent_mask | (word_type)1;
    static constexpr word_type inf_word = exponent_mask;
    static constexpr word_type zero_word = (word_type)0;
private:
    word_type word;
    struct make_from_bits_flag_t
    {
    };
    constexpr ieee754_soft_float(word_type word, make_from_bits_flag_t)
        : word(word)
    {
    }
    static constexpr int log_base_2(uintmax_t value)
    {
        return value == 0 ? -1 :
            (value >> 63) > 1U ? 64 + log_base_2((value >> 32) >> 32) :
            (value >> 32) > 0U ? 32 + log_base_2(value >> 32) :
            (value >> 16) > 0U ? 16 + log_base_2(value >> 16) :
            (value >> 8) > 0U ? 8 + log_base_2(value >> 8) :
            (value >> 4) > 0U ? 4 + log_base_2(value >> 4) :
            (value >> 2) > 0U ? 2 + log_base_2(value >> 2) :
            (value >> 1) > 0U ? 1 : 0;
    }
    static constexpr word_type construct_float_helper(word_type mantissa, word_type exponent)
    {
        return (mantissa & mantissa_mask) | ((exponent << exponent_shift) & exponent_mask);
    }
    static constexpr word_type construct_float_helper(std::pair<word_type, signed_type> v)
    {
        return std::get<0>(v) == (word_type)0 || std::get<1>(v) < min_denormal_normalized_exponent ? zero_word :
            std::get<1>(v) >= inf_nan_unbiased_exponent ? inf_word :
            std::get<1>(v) <= denormal_zero_unbiased_exponent ? construct_float_helper(std::get<0>(v) >> (intmax_t)(denormal_actual_exponent - std::get<1>(v)), denormal_zero_biased_exponent) :
            construct_float_helper(std::get<0>(v), (word_type)(std::get<1>(v) + (signed_type)exponent_bias));
    }
    static constexpr word_type construct_float_helper(bool is_negative, std::pair<word_type, signed_type> v)
    {
        return is_negative ? construct_float_helper(v) | sign_bit_mask : construct_float_helper(v);
    }
    static constexpr word_type construct_float_helper(std::pair<signed_type, signed_type> v)
    {
        return construct_float_helper(std::get<0>(v) < (signed_type)0, std::pair<word_type, signed_type>((word_type)abs(std::get<0>(v)), std::get<1>(v)));
    }
    template <size_t fractional_part_size = radix_point_shift, size_t N, typename = typename std::enable_if<(N >= bit_count)>::type>
    static constexpr std::pair<signed_type, signed_type> normalize(fixed_width_int<N> mantissa, signed_type exponent)
    {
        return mantissa == (signed_type)0 ? std::pair<signed_type, signed_type>((signed_type)mantissa, (signed_type)0) :
            (fixed_width_uint<N>)abs(mantissa) >> fractional_part_size == (fixed_width_uint<N>)0 ? normalize<fractional_part_size>(mantissa << 1, exponent - (signed_type)1) :
            (fixed_width_uint<N>)abs(mantissa) >> fractional_part_size > (fixed_width_uint<N>)1 ? normalize<fractional_part_size>(mantissa >> 1, exponent + (signed_type)1) :
            std::pair<signed_type, signed_type>((signed_type)mantissa, exponent);
    }
    template <size_t fractional_part_size = radix_point_shift, size_t N, typename = typename std::enable_if<(N < bit_count)>::type>
    static constexpr std::pair<word_type, signed_type> normalize(fixed_width_int<N> mantissa, signed_type exponent)
    {
        return normalize((signed_type)mantissa, exponent);
    }
    template <size_t fractional_part_size = radix_point_shift, size_t N, typename = typename std::enable_if<(N >= bit_count)>::type>
    static constexpr std::pair<word_type, signed_type> normalize(fixed_width_uint<N> mantissa, signed_type exponent)
    {
        return mantissa == (fixed_width_uint<N>)0 ? std::pair<word_type, signed_type>((word_type)mantissa, (signed_type)0) :
            mantissa >> fractional_part_size == (fixed_width_uint<N>)0 ? normalize<fractional_part_size>(mantissa << 1, exponent - (signed_type)1) :
            mantissa >> fractional_part_size > (fixed_width_uint<N>)1 ? normalize<fractional_part_size>(mantissa >> 1, exponent + (signed_type)1) :
            std::pair<word_type, signed_type>((word_type)mantissa, exponent);
    }
    template <size_t fractional_part_size = radix_point_shift, size_t N, typename = void, typename = typename std::enable_if<(N < bit_count)>::type>
    static constexpr std::pair<word_type, signed_type> normalize(fixed_width_uint<N> mantissa, signed_type exponent)
    {
        return normalize((word_type)mantissa, exponent);
    }
    template <size_t N>
    static constexpr word_type construct_float(fixed_width_int<N> mantissa, signed_type exponent)
    {
        return construct_float_helper(normalize(mantissa, exponent));
    }
    template <size_t N>
    static constexpr word_type construct_float(bool is_negative, fixed_width_uint<N> mantissa, signed_type exponent)
    {
        return construct_float_helper(is_negative, normalize(mantissa, exponent));
    }
    template <size_t N>
    static constexpr word_type int_to_float(fixed_width_int<N> value)
    {
        return construct_float(value, (signed_type)0);
    }
    template <size_t N>
    static constexpr word_type int_to_float(fixed_width_uint<N> value)
    {
        return construct_float(false, value, (signed_type)0);
    }
    constexpr word_type get_biased_exponent_value() const
    {
        return (word & exponent_mask) >> exponent_shift;
    }
    constexpr signed_type get_unbiased_exponent_value() const
    {
        return (signed_type)get_biased_exponent_value() - (signed_type)exponent_bias;
    }
    constexpr signed_type get_integer_exponent_value() const // exponent if the mantissa is interpreted as a integer
    {
        return (signed_type)get_biased_exponent_value() - ((signed_type)exponent_bias + (signed_type)radix_point_shift);
    }
    constexpr word_type get_mantissa_bits() const
    {
        return word & mantissa_mask;
    }
    constexpr word_type get_unsigned_mantissa_value() const
    {
        return get_unbiased_exponent_value() == denormal_zero_unbiased_exponent ? get_mantissa_bits() : get_mantissa_bits() | mantissa_implicit_bit;
    }
    constexpr signed_type get_signed_mantissa_value() const
    {
        return (word & sign_bit_mask) != (word_type)0 ? -(signed_type)get_unsigned_mantissa_value() : (signed_type)get_unsigned_mantissa_value();
    }
public:
    constexpr word_type get_bits() const
    {
        return word;
    }
    constexpr ieee754_soft_float()
        : word(0)
    {
    }
    static constexpr ieee754_soft_float from_bits(word_type word)
    {
        return ieee754_soft_float(word, make_from_bits_flag_t());
    }
    static constexpr ieee754_soft_float make_nan()
    {
        return from_bits(nan_word);
    }
    static constexpr ieee754_soft_float make_inf()
    {
        return from_bits(inf_word);
    }
    ieee754_soft_float(long double value)
    {
        word_type sign_bit = (word_type)0;
        if(std::signbit(value))
            sign_bit = sign_bit_mask;
        if(std::isnan(value))
        {
            word = sign_bit | nan_word;
        }
        else
        {
            if(sign_bit != (word_type)0)
                value = -value;
            if(std::isinf(value))
            {
                word = sign_bit | inf_word;
            }
            else if(value == 0)
            {
                word = sign_bit;
            }
            else if(value < std::pow(2.0L, 1 - (int)(uintmax_t)exponent_bias)) // denormal
            {
                int extra_shift = std::max<int>(0, (int)mantissa_bit_count - std::numeric_limits<long double>::digits);
                value = std::ldexp(value, (int)(uintmax_t)exponent_bias - 1 + mantissa_bit_count - extra_shift);
                word = ((word_type)std::floor(value) << extra_shift) | sign_bit;
            }
            else
            {
                int extra_shift = std::max<int>(0, (int)mantissa_bit_count - std::numeric_limits<long double>::digits);
                int exponent;
                value = std::frexp(value, &exponent);
                value = std::ldexp(value, (int)mantissa_bit_count - extra_shift);
                exponent--;
                if(exponent <= (int)(uintmax_t)exponent_bias)
                {
                    word = (((word_type)std::floor(value) << extra_shift) & mantissa_mask) | sign_bit | (((word_type)exponent + exponent_bias) << exponent_shift);
                }
                else
                {
                    word = sign_bit | inf_word;
                }
            }
        }
    }
    constexpr ieee754_soft_float(uintmax_t value)
        : word(int_to_float((word_type)value))
    {
    }
    constexpr ieee754_soft_float(intmax_t value)
        : word(int_to_float((signed_type)value))
    {
    }
    template <size_t N>
    constexpr ieee754_soft_float(fixed_width_uint<N> value)
        : word(int_to_float(value))
    {
    }
    template <size_t N>
    constexpr ieee754_soft_float(fixed_width_int<N> value)
        : word(int_to_float(value))
    {
    }
    constexpr bool isfinite() const
    {
        return (word & exponent_mask) != exponent_mask;
    }
    constexpr bool isnan() const
    {
        return !isfinite() && (word & mantissa_mask) != (word_type)0;
    }
    constexpr bool isinf() const
    {
        return !isfinite() && (word & mantissa_mask) == (word_type)0;
    }
    constexpr bool iszero() const
    {
        return get_biased_exponent_value() == denormal_zero_biased_exponent && get_mantissa_bits() == (word_type)0;
    }
    constexpr bool isdenormal() const
    {
        return get_biased_exponent_value() == denormal_zero_biased_exponent && get_mantissa_bits() != (word_type)0;
    }
    constexpr bool isnormal() const
    {
        return get_biased_exponent_value() != denormal_zero_biased_exponent && get_biased_exponent_value() != inf_nan_biased_exponent;
    }
    constexpr bool signbit() const
    {
        return (word & sign_bit_mask) != (word_type)0;
    }
private:
    template <typename T>
    static constexpr T constexpr_max(T a, T b)
    {
        return a < b ? b : a;
    }
    template <typename T>
    static constexpr T constexpr_min(T a, T b)
    {
        return a < b ? a : b;
    }
    template <size_t N>
    static constexpr fixed_width_int<N> round_to_zero_rshift(fixed_width_int<N> v, size_t shift_count)
    {
        return shift_count >= fixed_width_int<N>::bit_count ? (fixed_width_int<N>)0 : v < (fixed_width_int<N>)0 ? -(-v >> shift_count) : v >> shift_count;
    }
    static constexpr signed_type get_add_sub_compare_compute_at_exponent(bool a_is_negative, signed_type a_mantissa, signed_type a_exponent, bool b_is_negative, signed_type b_mantissa, signed_type b_exponent)
    {
        return a_mantissa == (signed_type)0 ? b_exponent :
            b_mantissa == (signed_type)0 ? a_mantissa :
                constexpr_max(a_exponent, b_exponent);
    }
    static constexpr signed_type get_add_mantissa(bool a_is_negative, signed_type a_mantissa, signed_type a_exponent, bool b_is_negative, signed_type b_mantissa, signed_type b_exponent)
    {
        return a_mantissa == (signed_type)0 ? b_mantissa :
            b_mantissa == (signed_type)0 ? a_mantissa :
            round_to_zero_rshift(a_mantissa, (intmax_t)(get_add_sub_compare_compute_at_exponent(a_is_negative, a_mantissa, a_exponent, b_is_negative, b_mantissa, b_exponent) - a_exponent)) + round_to_zero_rshift(b_mantissa, (intmax_t)(get_add_sub_compare_compute_at_exponent(a_is_negative, a_mantissa, a_exponent, b_is_negative, b_mantissa, b_exponent) - b_exponent));
    }
    static constexpr word_type add_helper(bool a_is_negative, signed_type a_mantissa, signed_type a_exponent, bool b_is_negative, signed_type b_mantissa, signed_type b_exponent)
    {
        return get_add_mantissa(a_is_negative, a_mantissa, a_exponent, b_is_negative, b_mantissa, b_exponent) == (signed_type)0 ? (a_is_negative && b_is_negative ? sign_bit_mask | zero_word : zero_word) : construct_float(get_add_mantissa(a_is_negative, a_mantissa, a_exponent, b_is_negative, b_mantissa, b_exponent), get_add_sub_compare_compute_at_exponent(a_is_negative, a_mantissa, a_exponent, b_is_negative, b_mantissa, b_exponent));
    }
public:
    friend constexpr ieee754_soft_float operator +(ieee754_soft_float a, ieee754_soft_float b)
    {
        return a.isnan() ? a :
            b.isnan() ? b :
            a.isinf() ?
                (b.isinf() ?
                    (a.signbit() == b.signbit() ? a :
                        (a.signbit() ? from_bits(nan_word | sign_bit_mask) : from_bits(nan_word | sign_bit_mask))) :
            a) :
            b.isinf() ? b :
            from_bits(add_helper(a.signbit(), a.get_signed_mantissa_value(), a.get_unbiased_exponent_value(), b.signbit(), b.get_signed_mantissa_value(), b.get_unbiased_exponent_value()));
    }
    friend constexpr ieee754_soft_float operator -(ieee754_soft_float a, ieee754_soft_float b)
    {
        return a + -b;
    }
    constexpr ieee754_soft_float operator -() const
    {
        return from_bits(word ^ sign_bit_mask);
    }
    constexpr ieee754_soft_float operator +() const
    {
        return *this;
    }
    friend constexpr ieee754_soft_float abs(ieee754_soft_float v)
    {
        return from_bits(v.word & ~sign_bit_mask);
    }
private:
    static constexpr word_type fractional_product(word_type a, word_type b)
    {
        return (word_type)((double_word_type)a * (double_word_type)b >> radix_point_shift);
    }
    static constexpr word_type fractional_quotient(word_type a, word_type b)
    {
        return (word_type)(((double_word_type)a << radix_point_shift) / (double_word_type)b);
    }
    static constexpr word_type product_helper(ieee754_soft_float a_abs, ieee754_soft_float b_abs)
    {
        return a_abs.isnan() ? a_abs.word :
            b_abs.isnan() ? b_abs.word :
            a_abs.isinf() ? (b_abs.iszero() ? nan_word : inf_word) :
            b_abs.isinf() ? (a_abs.iszero() ? nan_word : inf_word) :
            construct_float(false, fractional_product(a_abs.get_unsigned_mantissa_value(), b_abs.get_unsigned_mantissa_value()), a_abs.get_unbiased_exponent_value() + b_abs.get_unbiased_exponent_value());
    }
    static constexpr word_type quotient_helper(ieee754_soft_float a_abs, ieee754_soft_float b_abs)
    {
        return a_abs.isnan() ? a_abs.word :
            b_abs.isnan() ? b_abs.word :
            a_abs.isinf() ? (b_abs.isinf() ? nan_word : inf_word) :
            b_abs.isinf() ? zero_word :
            b_abs.iszero() ? (a_abs.iszero() ? nan_word : inf_word) :
            construct_float(false, fractional_quotient(a_abs.get_unsigned_mantissa_value(), b_abs.get_unsigned_mantissa_value()), a_abs.get_unbiased_exponent_value() - b_abs.get_unbiased_exponent_value());
    }
public:
    friend constexpr ieee754_soft_float operator *(ieee754_soft_float a, ieee754_soft_float b)
    {
        return from_bits(product_helper(abs(a), abs(b)) | ((a.word ^ b.word) & sign_bit_mask));
    }
    friend constexpr ieee754_soft_float operator /(ieee754_soft_float a, ieee754_soft_float b)
    {
        return from_bits(quotient_helper(abs(a), abs(b)) | ((a.word ^ b.word) & sign_bit_mask));
    }
    const ieee754_soft_float operator +=(ieee754_soft_float r)
    {
        return *this = *this + r;
    }
    const ieee754_soft_float operator -=(ieee754_soft_float r)
    {
        return *this = *this - r;
    }
    const ieee754_soft_float operator *=(ieee754_soft_float r)
    {
        return *this = *this * r;
    }
    const ieee754_soft_float operator /=(ieee754_soft_float r)
    {
        return *this = *this / r;
    }
    friend constexpr bool operator ==(ieee754_soft_float a, ieee754_soft_float b)
    {
        return a.isnan() ? false :
            a.iszero() ? b.iszero() :
            a.word == b.word;
    }
    friend constexpr bool operator !=(ieee754_soft_float a, ieee754_soft_float b)
    {
        return !(a == b);
    }
private:
    static constexpr int compare_helper(word_type a, word_type b)
    {
        return a == b ? 0 : a < b ? -1 : 1;
    }
    static constexpr int compare_helper(signed_type a, signed_type b)
    {
        return a == b ? 0 : a < b ? -1 : 1;
    }
    static constexpr int compare_helper(ieee754_soft_float a, ieee754_soft_float b)
    {
        return compare_helper(a.get_unsigned_mantissa_value() >> (uintmax_t)(constexpr_max(a.get_biased_exponent_value(), b.get_biased_exponent_value()) - a.get_biased_exponent_value()), b.get_unsigned_mantissa_value() >> (uintmax_t)(constexpr_max(a.get_biased_exponent_value(), b.get_biased_exponent_value()) - b.get_biased_exponent_value()));
    }
    static constexpr int compare(ieee754_soft_float a, ieee754_soft_float b, int unordered_value)
    {
        return a.isnan() || b.isnan() ? unordered_value :
            a.iszero() ? (b.iszero() ? 0 : (b.signbit() ? 1 : -1)) :
            b.iszero() || a.signbit() != b.signbit() ? (a.signbit() ? -1 : 1) :
            a.isinf() ? (b.isinf() ? 0 : (a.signbit() ? -1 : 1)) :
            a.signbit() ? -compare_helper(abs(a), abs(b)) : compare_helper(a, b);
    }
public:
    friend constexpr bool operator <(ieee754_soft_float a, ieee754_soft_float b)
    {
        return compare(a, b, 1) < 0;
    }
    friend constexpr bool operator >(ieee754_soft_float a, ieee754_soft_float b)
    {
        return b < a;
    }
    friend constexpr bool operator <=(ieee754_soft_float a, ieee754_soft_float b)
    {
        return compare(a, b, 1) <= 0;
    }
    friend constexpr bool operator >=(ieee754_soft_float a, ieee754_soft_float b)
    {
        return b <= a;
    }
private:
    template <size_t N, typename = typename std::enable_if<(N >= word_type::bit_count)>::type>
    constexpr fixed_width_int<N> convert_to_int() const
    {
        return get_integer_exponent_value() < (signed_type)0 ?
            round_to_zero_rshift((fixed_width_int<N>)get_signed_mantissa_value(), (intmax_t)(-get_integer_exponent_value())) :
            (fixed_width_int<N>)get_signed_mantissa_value() << (intmax_t)get_integer_exponent_value();
    }
    template <size_t N, typename = void, typename = typename std::enable_if<(N < word_type::bit_count)>::type>
    constexpr fixed_width_int<N> convert_to_int() const
    {
        return (fixed_width_int<N>)convert_to_int<word_type::bit_count>();
    }
    template <size_t N, typename = typename std::enable_if<(N >= word_type::bit_count)>::type>
    constexpr fixed_width_uint<N> convert_to_uint() const
    {
        return get_integer_exponent_value() < (signed_type)0 ?
            (fixed_width_uint<N>)get_unsigned_mantissa_value() >> (intmax_t)(-get_integer_exponent_value()) :
            (fixed_width_uint<N>)get_unsigned_mantissa_value() << (intmax_t)get_integer_exponent_value();
    }
    template <size_t N, typename = void, typename = typename std::enable_if<(N < word_type::bit_count)>::type>
    constexpr fixed_width_uint<N> convert_to_uint() const
    {
        return (fixed_width_uint<N>)convert_to_uint<word_type::bit_count>();
    }
public:
    template <size_t N>
    explicit constexpr operator fixed_width_int<N>() const
    {
        return isnan() ? (fixed_width_int<N>)0 :
            isinf() ? (signbit() ? std::numeric_limits<fixed_width_int<N>>::min() : std::numeric_limits<fixed_width_int<N>>::max()) :
            *this >= -(ieee754_soft_float)std::numeric_limits<fixed_width_int<N>>::min() ? std::numeric_limits<fixed_width_int<N>>::max() :
            *this <= (ieee754_soft_float)std::numeric_limits<fixed_width_int<N>>::min() ? std::numeric_limits<fixed_width_int<N>>::min() :
            *this < (ieee754_soft_float)1 && *this > (ieee754_soft_float)-1 ? (fixed_width_int<N>)0 :
            convert_to_int<N>();
    }
    template <size_t N>
    explicit constexpr operator fixed_width_uint<N>() const
    {
        return isnan() ? (fixed_width_uint<N>)0 :
            isinf() ? (signbit() ? (fixed_width_uint<N>)0 : std::numeric_limits<fixed_width_uint<N>>::max()) :
            *this >= (ieee754_soft_float)((fixed_width_uint<2 * N>)1 << N) ? std::numeric_limits<fixed_width_uint<N>>::max() :
            *this < (ieee754_soft_float)1 ? (fixed_width_uint<N>)0 :
            convert_to_uint<N>();
    }
    friend int main(int argc, char **argv);
    #warning remove friend main
};

typedef ieee754_soft_float<24, 8, true> ieee754_soft_float_32;
typedef ieee754_soft_float<53, 11, true> ieee754_soft_float_64;
typedef ieee754_soft_float<64, 15, false, 128 - 80> ieee754_soft_float_80;
typedef ieee754_soft_float<113, 15, true> ieee754_soft_float_128;
typedef ieee754_soft_float<237, 19, true> ieee754_soft_float_256;
typedef ieee754_soft_float<489, 23, true> ieee754_soft_float_512;
typedef ieee754_soft_float<997, 27, true> ieee754_soft_float_1024;
typedef ieee754_soft_float<2048 - 31 - 1 + 1, 31, true> ieee754_soft_float_2048;
typedef ieee754_soft_float<4096 - 35 - 1 + 1, 35, true> ieee754_soft_float_4096;
typedef ieee754_soft_float<8192 - 39 - 1 + 1, 39, true> ieee754_soft_float_8192;

template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::sign_bit_mask;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::exponent_mask;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::exponent_bias;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::inf_nan_biased_exponent;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::denormal_zero_biased_exponent;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::signed_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::inf_nan_unbiased_exponent;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::signed_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::denormal_zero_unbiased_exponent;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::signed_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::denormal_actual_exponent;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::signed_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::min_denormal_normalized_exponent;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::mantissa_mask;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::mantissa_implicit_bit;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::nan_word;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::inf_word;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::zero_word;


#endif // FLOAT_TYPES_H_INCLUDED
