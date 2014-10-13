#ifndef FLOAT_TYPES_H_INCLUDED
#define FLOAT_TYPES_H_INCLUDED

#include "int_types.h"
#include <cmath>
#include <type_traits>

#if __cplusplus > 201103L
#define cpp14_constexpr constexpr
#else
#define cpp14_constexpr
#endif

template <size_t N, typename = void>
struct ieee754_soft_float_std;

template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize = 0>
class ieee754_soft_float
{
    template <size_t MBC, size_t EBC, bool MMII, size_t PS>
    friend class ieee754_soft_float;
public:
    static constexpr size_t mantissa_bit_count = mantissaBitCount;
    static_assert(mantissa_bit_count >= 2, "mantissa too small");
    static constexpr size_t exponent_bit_count = exponentBitCount;
    static_assert(mantissa_bit_count >= 2, "exponent too small");
    static constexpr bool mantissa_msb_is_implicit = mantissaMsbIsImplicit;
    static constexpr size_t padding_size = paddingSize;
    static constexpr size_t bit_count = mantissa_bit_count + exponent_bit_count + 1 - (mantissa_msb_is_implicit ? 1 : 0);
    static constexpr size_t padded_bit_count = bit_count + padding_size;
    typedef fixed_width_uint<padded_bit_count> word_type;
    typedef fixed_width_int<padded_bit_count> signed_type;
    typedef fixed_width_uint<padded_bit_count * 2> double_word_type;
private:
    static constexpr size_t sign_bit_shift = bit_count - 1;
    static constexpr size_t exponent_shift = sign_bit_shift - exponent_bit_count;
    static constexpr size_t mantissa_shift = 0;
    static constexpr size_t radix_point_shift = mantissa_bit_count - 1;
    static constexpr word_type sign_bit_mask = (word_type)1 << sign_bit_shift;
    static constexpr word_type exponent_mask = (((word_type)1 << exponent_bit_count) - (word_type)1) << exponent_shift;
    static constexpr word_type exponent_bias = ((word_type)1 << (exponent_bit_count - 1)) - (word_type)1;
    static constexpr word_type inf_nan_exponent = exponent_mask >> exponent_shift;
    static constexpr word_type inf_nan_biased_exponent = inf_nan_exponent;
    static constexpr word_type denormal_zero_exponent = (word_type)0;
    static constexpr signed_type inf_nan_unbiased_exponent = (signed_type)inf_nan_biased_exponent - (signed_type)exponent_bias;
    static constexpr word_type denormal_biased_actual_exponent = (word_type)1;
    static constexpr signed_type denormal_unbiased_actual_exponent = (signed_type)denormal_biased_actual_exponent - (signed_type)exponent_bias;
    static constexpr signed_type min_denormal_normalized_exponent = denormal_unbiased_actual_exponent - (signed_type)(intmax_t)mantissa_bit_count;
    static constexpr word_type mantissa_mask = ((word_type)1 << (mantissa_bit_count - (mantissa_msb_is_implicit ? 1 : 0))) - (word_type)1;
    static constexpr word_type mantissa_implicit_bit = mantissa_msb_is_implicit ? (word_type)1 << (mantissa_bit_count - 1) : (word_type)0;
    static constexpr word_type nan_word = exponent_mask | (word_type)1;
    static constexpr word_type inf_word = exponent_mask;
    static constexpr word_type zero_word = (word_type)0;
    static constexpr word_type one_word = ((word_type)1 << (mantissa_bit_count - 1)) | (exponent_bias << exponent_shift);
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
            std::get<1>(v) < denormal_unbiased_actual_exponent ? construct_float_helper(std::get<0>(v) >> (intmax_t)(denormal_unbiased_actual_exponent - std::get<1>(v)), denormal_zero_exponent) :
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
    static constexpr std::pair<signed_type, signed_type> normalize_signed_helper(bool is_negative, std::pair<word_type, signed_type> v)
    {
        return is_negative ? std::pair<signed_type, signed_type>(-(signed_type)std::get<0>(v), std::get<1>(v)) : std::pair<signed_type, signed_type>((signed_type)std::get<0>(v), std::get<1>(v));
    }
    template <size_t fractional_part_size = radix_point_shift, size_t N>
    static constexpr std::pair<signed_type, signed_type> normalize(fixed_width_int<N> mantissa, signed_type exponent)
    {
        return normalize_signed_helper(mantissa < (fixed_width_int<N>)0, normalize<fractional_part_size>((fixed_width_uint<N>)mantissa, exponent));
    }
    template <size_t fractional_part_size = radix_point_shift, size_t N, typename = typename std::enable_if<(N >= bit_count)>::type>
    static constexpr std::pair<word_type, signed_type> normalize(fixed_width_uint<N> mantissa, signed_type exponent)
    {
#if __cplusplus > 201103L
        if(mantissa == (fixed_width_uint<N>)0)
            return std::pair<word_type, signed_type>((word_type)mantissa, (signed_type)0);
        intmax_t lshift_amount = (intmax_t)fractional_part_size - ilog2(mantissa);
        exponent -= (signed_type)lshift_amount;
        if(lshift_amount > 0)
            mantissa <<= lshift_amount;
        else if(lshift_amount < 0)
            mantissa >>= -lshift_amount;
        return std::pair<word_type, signed_type>((word_type)mantissa, exponent);
#else
        return mantissa == (fixed_width_uint<N>)0 ? std::pair<word_type, signed_type>((word_type)mantissa, (signed_type)0) :
            mantissa >> fractional_part_size == (fixed_width_uint<N>)0 ? normalize<fractional_part_size>(mantissa << 1, exponent - (signed_type)1) :
            mantissa >> fractional_part_size > (fixed_width_uint<N>)1 ? normalize<fractional_part_size>(mantissa >> 1, exponent + (signed_type)1) :
            std::pair<word_type, signed_type>((word_type)mantissa, exponent);
#endif
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
        return construct_float(value, (signed_type)radix_point_shift);
    }
    template <size_t N>
    static constexpr word_type int_to_float(fixed_width_uint<N> value)
    {
        return construct_float(false, value, (signed_type)radix_point_shift);
    }
    constexpr word_type get_exponent_bits() const
    {
        return (word & exponent_mask) >> exponent_shift;
    }
    constexpr word_type get_biased_exponent_value() const
    {
        return get_exponent_bits() == denormal_zero_exponent ? denormal_biased_actual_exponent : get_exponent_bits();
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
        return get_exponent_bits() == denormal_zero_exponent ? get_mantissa_bits() : get_mantissa_bits() | mantissa_implicit_bit;
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
    static constexpr ieee754_soft_float make_denormal_min()
    {
        return from_bits((word_type)1);
    }
    static constexpr ieee754_soft_float make_min()
    {
        return from_bits((word_type)1 << exponent_shift);
    }
    static constexpr ieee754_soft_float make_max()
    {
        return from_bits(((exponent_mask - (word_type)1) & exponent_mask) | mantissa_mask);
    }
    static constexpr ieee754_soft_float make_epsilon()
    {
        return from_bits(construct_float((signed_type)1 << radix_point_shift, -(signed_type)mantissa_bit_count));
    }
    static constexpr ieee754_soft_float make_round_error()
    {
        return from_bits(construct_float((signed_type)1 << radix_point_shift, (signed_type)-1));
    }
    explicit ieee754_soft_float(long double value)
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
            else if(exponent_bias < (word_type)0x1000000 && value < std::pow(2.0L, 1 - (int)(uintmax_t)exponent_bias)) // denormal
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
                if((signed_type)(intmax_t)exponent <= (signed_type)exponent_bias)
                {
                    word = (((word_type)(uintmax_t)std::floor(value) << extra_shift) & mantissa_mask) | sign_bit | (((word_type)exponent + exponent_bias) << exponent_shift);
                }
                else
                {
                    word = sign_bit | inf_word;
                }
            }
        }
    }
    explicit constexpr ieee754_soft_float(uintmax_t value)
        : word(int_to_float((word_type)value))
    {
    }
    explicit constexpr ieee754_soft_float(intmax_t value)
        : word(int_to_float((signed_type)value))
    {
    }
    template <size_t N>
    explicit constexpr ieee754_soft_float(fixed_width_uint<N> value)
        : word(int_to_float(value))
    {
    }
    template <size_t N>
    explicit constexpr ieee754_soft_float(fixed_width_int<N> value)
        : word(int_to_float(value))
    {
    }
private:
    template <size_t r_mantissa_bit_count, size_t N, typename = typename std::enable_if<(r_mantissa_bit_count > mantissa_bit_count)>::type>
    static constexpr word_type convert_mantissa(fixed_width_uint<N> v)
    {
        return (word_type)(v >> (r_mantissa_bit_count - mantissa_bit_count));
    }
    template <size_t r_mantissa_bit_count, size_t N, typename = void, typename = typename std::enable_if<(r_mantissa_bit_count < mantissa_bit_count)>::type>
    static constexpr word_type convert_mantissa(fixed_width_uint<N> v)
    {
        return (word_type)v << (mantissa_bit_count - r_mantissa_bit_count);
    }
    template <size_t r_mantissa_bit_count, size_t N, typename = void, typename = void, typename = typename std::enable_if<(r_mantissa_bit_count == mantissa_bit_count)>::type>
    static constexpr word_type convert_mantissa(fixed_width_uint<N> v)
    {
        return (word_type)v;
    }
    static constexpr word_type convert_nan_mantissa_helper(word_type v)
    {
        return v == (word_type)0 ? (word_type)1 : v;
    }
    template <size_t r_mantissa_bit_count, size_t N>
    static constexpr word_type convert_nan_mantissa(fixed_width_uint<N> v)
    {
        return convert_nan_mantissa_helper(convert_mantissa<r_mantissa_bit_count>(v));
    }
    template <size_t N, typename = typename std::enable_if<(N > signed_type::bit_count)>::type>
    static constexpr bool exponent_too_big(fixed_width_int<N> v)
    {
        return v >= (fixed_width_int<N>)1 << exponent_bit_count;
    }
    template <size_t N, typename = void, typename = typename std::enable_if<(N <= signed_type::bit_count)>::type>
    static constexpr bool exponent_too_big(fixed_width_int<N> v)
    {
        return (signed_type)v >= (signed_type)1 << exponent_bit_count;
    }
    template <size_t N, typename = typename std::enable_if<(N > signed_type::bit_count)>::type>
    static constexpr bool exponent_too_small(fixed_width_int<N> v)
    {
        return v < (fixed_width_int<N>)-1 << exponent_bit_count;
    }
    template <size_t N, typename = void, typename = typename std::enable_if<(N <= signed_type::bit_count)>::type>
    static constexpr bool exponent_too_small(fixed_width_int<N> v)
    {
        return (signed_type)v < (signed_type)-1 << exponent_bit_count;
    }
    template <typename r_float_type>
    static constexpr word_type convert_float_no_sign(r_float_type v)
    {
        return v.isnan() ? exponent_mask | (convert_nan_mantissa<r_float_type::mantissa_bit_count>(v.get_unsigned_mantissa_value()) & mantissa_mask) :
            v.isinf() ? inf_word :
            v.iszero() || exponent_too_small(v.get_unbiased_exponent_value()) ? zero_word :
            exponent_too_big(v.get_unbiased_exponent_value()) ? inf_word :
            construct_float(false, v.get_unsigned_mantissa_value(), (signed_type)v.get_unbiased_exponent_value() + (signed_type)((intmax_t)radix_point_shift - (intmax_t)r_float_type::radix_point_shift));
    }
public:
    template <size_t MBC, size_t EBC, bool MMII, size_t PS>
    explicit constexpr ieee754_soft_float(ieee754_soft_float<MBC, EBC, MMII, PS> value)
        : word(value.signbit() ? sign_bit_mask | convert_float_no_sign(value) : convert_float_no_sign(value))
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
        return get_exponent_bits() == denormal_zero_exponent && get_mantissa_bits() == (word_type)0;
    }
    constexpr bool isdenormal() const
    {
        return get_exponent_bits() == denormal_zero_exponent && get_mantissa_bits() != (word_type)0;
    }
    constexpr bool isnormal() const
    {
        return get_exponent_bits() != denormal_zero_exponent && get_exponent_bits() != inf_nan_exponent;
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
    static constexpr fixed_width_int<N> round_to_zero_rshift(fixed_width_int<N> v, intmax_t shift_count)
    {
        return shift_count >= fixed_width_int<N>::bit_count ? (fixed_width_int<N>)0 :
            shift_count < 0 ? v << -shift_count :
            v < (fixed_width_int<N>)0 ? -(-v >> shift_count) : v >> shift_count;
    }
    template <size_t N>
    static constexpr fixed_width_uint<N> round_to_zero_rshift(fixed_width_uint<N> v, intmax_t shift_count)
    {
        return shift_count >= fixed_width_int<N>::bit_count ? (fixed_width_uint<N>)0 :
            shift_count < 0 ? v << -shift_count :
            v >> shift_count;
    }
    static constexpr signed_type get_add_sub_compare_compute_at_exponent(bool a_is_negative, signed_type a_mantissa, signed_type a_exponent, bool b_is_negative, signed_type b_mantissa, signed_type b_exponent)
    {
        return a_mantissa == (signed_type)0 ? b_exponent :
            b_mantissa == (signed_type)0 ? a_mantissa :
                constexpr_max(a_exponent, b_exponent) - (signed_type)2;
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
    static constexpr word_type round_word(word_type v)
    {
        return (v & (word_type)1) + (v >> 1);
    }
    static constexpr word_type fractional_quotient(word_type a, word_type b)
    {
        return round_word((word_type)(((double_word_type)a << (radix_point_shift + 1)) / (double_word_type)b));
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
    cpp14_constexpr const ieee754_soft_float operator +=(ieee754_soft_float r)
    {
        return *this = *this + r;
    }
    cpp14_constexpr const ieee754_soft_float operator -=(ieee754_soft_float r)
    {
        return *this = *this - r;
    }
    cpp14_constexpr const ieee754_soft_float operator *=(ieee754_soft_float r)
    {
        return *this = *this * r;
    }
    cpp14_constexpr const ieee754_soft_float operator /=(ieee754_soft_float r)
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
            *this < from_bits(one_word) && *this > -from_bits(one_word) ? (fixed_width_int<N>)0 :
            convert_to_int<N>();
    }
    template <size_t N>
    explicit constexpr operator fixed_width_uint<N>() const
    {
        return isnan() ? (fixed_width_uint<N>)0 :
            isinf() ? (signbit() ? (fixed_width_uint<N>)0 : std::numeric_limits<fixed_width_uint<N>>::max()) :
            *this >= (ieee754_soft_float)((fixed_width_uint<2 * N>)1 << N) ? std::numeric_limits<fixed_width_uint<N>>::max() :
            *this < (ieee754_soft_float)(word_type)1 ? (fixed_width_uint<N>)0 :
            convert_to_uint<N>();
    }
private:
    constexpr word_type get_floor_mask() const
    {
        return ((mantissa_mask | mantissa_implicit_bit) << (intmax_t)-get_integer_exponent_value()) & (mantissa_mask | mantissa_implicit_bit);
    }
    constexpr word_type floor_helper() const
    {
        return get_integer_exponent_value() >= (signed_type)0 ? get_bits() :
                get_integer_exponent_value() < -(signed_type)mantissa_bit_count ? (signbit() ? one_word | sign_bit_mask : zero_word) :
                (((mantissa_mask | mantissa_implicit_bit) ^ get_floor_mask()) & get_unsigned_mantissa_value()) == (word_type)0 ? get_bits() :
                signbit() ? (from_bits(construct_float(true, get_unsigned_mantissa_value() & get_floor_mask(), get_unbiased_exponent_value())) - from_bits(one_word)).get_bits() : construct_float(false, get_unsigned_mantissa_value() & get_floor_mask(), get_unbiased_exponent_value());
    }
public:
    friend constexpr ieee754_soft_float floor(ieee754_soft_float v)
    {
        return v.isnan() || v.isinf() || v.iszero() ? v :
            v < from_bits(one_word) && v >= from_bits(zero_word) ? from_bits(zero_word) : from_bits(v.floor_helper());
    }
    friend constexpr ieee754_soft_float ceil(ieee754_soft_float v)
    {
        return -floor(-v);
    }
    friend constexpr ieee754_soft_float trunc(ieee754_soft_float v)
    {
        return v.signbit() ? ceil(v) : floor(v);
    }
private:
    static constexpr word_type fractional_sqrt(word_type v)
    {
        return round_word((word_type)floor_sqrt((double_word_type)v << (radix_point_shift + 2)));
    }
public:
    friend constexpr ieee754_soft_float sqrt(ieee754_soft_float v)
    {
        return v.isnan() || v.iszero() ? v :
            v.signbit() ? from_bits(nan_word) :
            v.isinf() ? from_bits(inf_word) :
            (v.get_unbiased_exponent_value() & (signed_type)1) != (signed_type)0 ?
                from_bits(construct_float(false, fractional_sqrt(v.get_unsigned_mantissa_value() << 1), (v.get_unbiased_exponent_value() - (signed_type)1) / (signed_type)2)) :
            from_bits(construct_float(false, fractional_sqrt(v.get_unsigned_mantissa_value()), v.get_unbiased_exponent_value() / (signed_type)2));
    }
    static ieee754_soft_float sqrt_2()
    {
        static ieee754_soft_float retval = sqrt((ieee754_soft_float)(word_type)2);
        return retval;
    }
    static ieee754_soft_float sqrt_1_2()
    {
        static ieee754_soft_float retval = sqrt((ieee754_soft_float)(word_type)1 / (ieee754_soft_float)(word_type)2);
        return retval;
    }
private:
    typedef typename ieee754_soft_float_std<padded_bit_count < 16 ? 32 : padded_bit_count * 2>::type double_precision_type;
    static constexpr double_precision_type square(double_precision_type v)
    {
        return v * v;
    }
    static ieee754_soft_float calc_pi_helper(double_precision_type a, double_precision_type b, double_precision_type t, double_precision_type p)
    {
        return a == b ? (ieee754_soft_float)(a * a / t) :
            calc_pi_helper((a + b) / (double_precision_type)(fixed_width_uint<padded_bit_count * 2>)2, sqrt(a * b), t - p * square(a - (a + b) / (double_precision_type)(fixed_width_uint<padded_bit_count * 2>)2), p * (double_precision_type)(fixed_width_uint<padded_bit_count * 2>)2);
    }
    static ieee754_soft_float calc_pi()
    {
        return calc_pi_helper((double_precision_type)(fixed_width_uint<padded_bit_count * 2>)1, double_precision_type::sqrt_1_2(), (double_precision_type)(fixed_width_uint<padded_bit_count * 2>)1 / (double_precision_type)(fixed_width_uint<padded_bit_count * 2>)4, (double_precision_type)(fixed_width_uint<padded_bit_count * 2>)1);
    }
public:
    static ieee754_soft_float pi()
    {
        static ieee754_soft_float retval = calc_pi();
        return retval;
    };
private:
    static double_precision_type calc_log1p_range_limited(double_precision_type v)
    {
        if(abs(v) >= (double_precision_type)(uintmax_t)1 / (double_precision_type)(uintmax_t)256)
        {
            return (double_precision_type)(uintmax_t)2 * calc_log1p_range_limited(sqrt(v + (double_precision_type)(uintmax_t)1) - (double_precision_type)(uintmax_t)1);
        }
        double_precision_type power = v;
        double_precision_type retval = (double_precision_type)(uintmax_t)0;
        for(uintmax_t i = 1; ; i++, power *= v)
        {
            double_precision_type divisor(i);
            double_precision_type term = power / divisor;
            double_precision_type old_retval = retval;
            if(i % 2 == 1)
                retval += term;
            else
                retval -= term;
            if(retval == old_retval)
                break;
        }
        return retval;
    }
    static double_precision_type ln2_double()
    {
        static double_precision_type retval = calc_log1p_range_limited((double_precision_type)(uintmax_t)1);
        return retval;
    }
    static cpp14_constexpr ieee754_soft_float log2_helper(ieee754_soft_float v)
    {
        signed_type exponent = ilogb(v);
        ieee754_soft_float retval = (ieee754_soft_float)exponent;
        v = scalbn(v, -exponent);
        for(size_t i = 0; i < mantissa_bit_count; i++)
        {
            v *= v;
            if(v >= (ieee754_soft_float)(uintmax_t)2)
            {
                retval += scalbn((ieee754_soft_float)(uintmax_t)1, (signed_type)(-1 - (intmax_t)i));
                v = scalbn(v, (signed_type)(intmax_t)-1);
            }
        }
        return retval;
    }
public:
    static constexpr signed_type ilogb_zero = std::numeric_limits<signed_type>::min();
    static constexpr signed_type ilogb_nan = -std::numeric_limits<signed_type>::max();
    friend constexpr signed_type ilogb(ieee754_soft_float v)
    {
        return v.isnan() ? ilogb_nan :
            v.iszero() ? ilogb_zero :
            v.isinf() ? std::numeric_limits<signed_type>::max() :
            std::get<1>(normalize(v.get_unsigned_mantissa_value(), v.get_unbiased_exponent_value()));
    }
    friend constexpr ieee754_soft_float scalbn(ieee754_soft_float v, signed_type exponent)
    {
        return v.isnan() || v.isinf() || v.iszero() ? v :
            exponent > (signed_type)exponent_bias * (signed_type)(intmax_t)4 ? (v.signbit() ? -make_inf() : make_inf()) :
            exponent < (signed_type)exponent_bias * (signed_type)(intmax_t)-4 ? (v.signbit() ? -from_bits(zero_word) : from_bits(zero_word)) :
            from_bits(construct_float(v.signbit(), v.get_unsigned_mantissa_value(), v.get_unbiased_exponent_value() + exponent));
    }
    friend constexpr ieee754_soft_float ldexp(ieee754_soft_float v, signed_type exponent)
    {
        return scalbn(v, exponent);
    }
    friend cpp14_constexpr ieee754_soft_float log2(ieee754_soft_float v)
    {
        if(v.isnan())
            return v;
        if(v.signbit() || v.iszero())
            return make_nan();
        if(v.isinf())
            return v;
        auto exponent = ilogb(v);
        return (ieee754_soft_float)exponent + log2_helper(scalbn(v, -exponent));
    }
    static cpp14_constexpr ieee754_soft_float log2_10()
    {
#if __cplusplus > 201103L
        constexpr
#else
        static
#endif
        ieee754_soft_float retval = (ieee754_soft_float)log2((double_precision_type)(uintmax_t)10);
        return retval;
    }
    static cpp14_constexpr ieee754_soft_float log10_2()
    {
#if __cplusplus > 201103L
        constexpr
#else
        static
#endif
        ieee754_soft_float retval = (ieee754_soft_float)((double_precision_type)(uintmax_t)1 / log2((double_precision_type)(uintmax_t)10));
        return retval;
    }
    friend cpp14_constexpr ieee754_soft_float log10(ieee754_soft_float v)
    {
        return log2(v) / log10_2();
    }
    static ieee754_soft_float ln2()
    {
        return (ieee754_soft_float)ln2_double();
    }
    friend int main(int argc, char **argv);
    #warning remove friend main
};

template <size_t N>
struct ieee754_soft_float_std<N, typename std::enable_if<N >= 64 && (N & (N - 1)) == 0>::type>
{
private:
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
    static constexpr size_t exponent_bit_count = 4 * log_base_2(N) - 13;
public:
    typedef ieee754_soft_float<N - exponent_bit_count + 1 - 1, exponent_bit_count, true> type;
};

template <>
struct ieee754_soft_float_std<8, void>
{
    typedef ieee754_soft_float<4, 4, true> type;
};

template <>
struct ieee754_soft_float_std<16, void>
{
    typedef ieee754_soft_float<11, 5, true> type;
};

template <>
struct ieee754_soft_float_std<32, void>
{
    typedef ieee754_soft_float<24, 8, true> type;
};

typedef typename ieee754_soft_float_std<8>::type ieee754_soft_float_8;
typedef typename ieee754_soft_float_std<16>::type ieee754_soft_float_16;
typedef typename ieee754_soft_float_std<32>::type ieee754_soft_float_32;
typedef typename ieee754_soft_float_std<64>::type ieee754_soft_float_64;
typedef typename ieee754_soft_float_std<128>::type ieee754_soft_float_128;
typedef typename ieee754_soft_float_std<256>::type ieee754_soft_float_256;
typedef typename ieee754_soft_float_std<512>::type ieee754_soft_float_512;
typedef typename ieee754_soft_float_std<1024>::type ieee754_soft_float_1k;
typedef typename ieee754_soft_float_std<2048>::type ieee754_soft_float_2k;
typedef typename ieee754_soft_float_std<4096>::type ieee754_soft_float_4k;
typedef typename ieee754_soft_float_std<8192>::type ieee754_soft_float_8k;
typedef typename ieee754_soft_float_std<16384>::type ieee754_soft_float_16k;
typedef typename ieee754_soft_float_std<32768>::type ieee754_soft_float_32k;
typedef ieee754_soft_float<64, 15, false, 128 - 80> ieee754_soft_float_i387_80;

template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::sign_bit_mask;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::exponent_mask;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::exponent_bias;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::inf_nan_biased_exponent;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::denormal_biased_actual_exponent;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::denormal_zero_exponent;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::signed_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::inf_nan_unbiased_exponent;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::signed_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::denormal_unbiased_actual_exponent;
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
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::word_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::one_word;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::signed_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::ilogb_nan;
template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
constexpr typename ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::signed_type ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::ilogb_zero;
//template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
//constexpr ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize> ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::pi;
//template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
//constexpr ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize> ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::sqrt_2;
//template <size_t mantissaBitCount, size_t exponentBitCount, bool mantissaMsbIsImplicit, size_t paddingSize>
//constexpr ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize> ieee754_soft_float<mantissaBitCount, exponentBitCount, mantissaMsbIsImplicit, paddingSize>::sqrt_1_2;

namespace std
{
template <size_t MBC, size_t EBC, bool MMII, size_t PS>
struct numeric_limits<ieee754_soft_float<MBC, EBC, MMII, PS>>
{
private:
    typedef ieee754_soft_float<MBC, EBC, MMII, PS> fp_type;
private:
    static constexpr intmax_t constexpr_ifloor(long double v)
    {
        return v < (intmax_t)v ? (intmax_t)v - 1 : (intmax_t)v;
    }
    static constexpr intmax_t constexpr_iceil(long double v)
    {
        return v > (intmax_t)v ? (intmax_t)v + 1 : (intmax_t)v;
    }
    static constexpr long double log10_2 = 0.30102999566398119521373889472449302676818988L;
public:
    static constexpr bool is_specialized = true;
    static constexpr fp_type min()
    {
        return fp_type::make_min();
    }
    static constexpr fp_type max()
    {
        return fp_type::make_max();
    }
    static constexpr fp_type lowest()
    {
        return -fp_type::make_max();
    }
    static constexpr int digits = fp_type::mantissa_bit_count;
    static constexpr int digits10 = (int)constexpr_ifloor(log10_2 * digits);
    static constexpr int max_digits10 = (int)constexpr_iceil(log10_2 * digits);
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = false;
    static constexpr bool is_exact = false;
    static constexpr int radix = 2;
    static constexpr fp_type epsilon()
    {
        return fp_type::make_epsilon();
    }
    static constexpr fp_type round_error()
    {
        return fp_type::make_round_error();
    }
    static constexpr intmax_t max_exponent = (intmax_t)1 << (fp_type::exponent_bit_count - 1);
    static constexpr intmax_t min_exponent = (intmax_t)1 - max_exponent;
    static constexpr intmax_t min_exponent10 = constexpr_iceil(min_exponent * log10_2);
    static constexpr intmax_t max_exponent10 = constexpr_ifloor(max_exponent * log10_2);
    static constexpr bool has_infinity = true;
    static constexpr bool has_quiet_NaN = true;
    static constexpr bool has_signaling_NaN = true;
    static constexpr float_denorm_style has_denorm = float_denorm_style::denorm_present;
    static constexpr bool has_denorm_loss = false;
    static constexpr fp_type infinity()
    {
        return fp_type::make_inf();
    }
    static constexpr fp_type quiet_NaN()
    {
        return fp_type::make_nan();
    }
    static constexpr fp_type signaling_NaN()
    {
        return fp_type::make_nan();
    }
    static constexpr fp_type denorm_min()
    {
        return fp_type::make_denormal_min();
    }
    static constexpr bool is_iec559 = true;
    static constexpr bool is_bounded = false;
    static constexpr bool is_modulo = false;
    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;
    static constexpr float_round_style round_style = float_round_style::round_to_nearest;
};
}

#endif // FLOAT_TYPES_H_INCLUDED
