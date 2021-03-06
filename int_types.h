#ifndef INT_TYPES_H_INCLUDED
#define INT_TYPES_H_INCLUDED

#include <cstdint>
#include <limits>
#include <cassert>
#include <tuple>
#include <istream>
#include <ostream>
#include <string>
#include <sstream>
#include <type_traits>
#include <stdexcept>

#if __cplusplus > 201103L
#define cpp14_constexpr constexpr
#else
#define cpp14_constexpr inline
#endif

template <size_t bitCount, typename = void>
struct fixed_width_uint_multiply_helper;

template <size_t bitCount, typename = void>
struct fixed_width_int_multiply_helper;

template <size_t bitCount>
struct fixed_width_uint final
{
    template <size_t N>
    friend struct fixed_width_uint;
    static constexpr size_t bit_count = bitCount;
    static_assert(bit_count >= 8, "invalid bit count for fixed_width_uint : must be at least 8");
    static_assert((bit_count & (bit_count - 1)) == 0, "invalid bit count for fixed_width_uint : must be a power of 2");
private:
    typedef fixed_width_uint<bit_count / 2> half_width_type;
    typedef fixed_width_uint<bit_count / 4> quarter_width_type;
    static_assert(std::numeric_limits<std::uintmax_t>::radix == 2, "radix of uintmax_t is not 2");
    static constexpr size_t uintmax_bit_count = (size_t)std::numeric_limits<std::uintmax_t>::digits;
    half_width_type lowPart;
    half_width_type highPart;
    static constexpr half_width_type quarter_mask = ((half_width_type)1 << quarter_width_type::bit_count) - (half_width_type)1;
    static constexpr bool get_sum_carry(half_width_type a, half_width_type b)
    {
        return (a + b) < a;
    }
    static constexpr bool get_difference_borrow(half_width_type a, half_width_type b)
    {
        return a < b;
    }
    constexpr fixed_width_uint(half_width_type highPart, half_width_type lowPart)
        : lowPart(lowPart), highPart(highPart)
    {
    }
    static constexpr int ilog_base_2(fixed_width_uint v, size_t current_bit_count = bit_count / 2)
    {
        return v == (fixed_width_uint)0 ? -1 :
            current_bit_count >= 1 ? ((v >> current_bit_count) > (fixed_width_uint)0 ? current_bit_count + ilog_base_2(v >> current_bit_count, current_bit_count / 2) : ilog_base_2(v, current_bit_count / 2)) :
            0;
    }
    static constexpr fixed_width_uint get_double_product(half_width_type a, half_width_type b)
    {
        return wide_multiply(a, b);
    }
    static constexpr fixed_width_uint get_product(fixed_width_uint a, fixed_width_uint b)
    {
        return a.highPart == (half_width_type)0 && b.highPart == (half_width_type)0 ? get_double_product(a.lowPart, b.lowPart) :
            get_double_product(a.lowPart, b.lowPart) + fixed_width_uint(a.lowPart * b.highPart + a.highPart * b.lowPart, (half_width_type)0);
    }
    static constexpr fixed_width_uint get_sum(fixed_width_uint a, fixed_width_uint b)
    {
        return fixed_width_uint(get_sum_carry(a.lowPart, b.lowPart) ? a.highPart + b.highPart + (half_width_type)1 : a.highPart + b.highPart, a.lowPart + b.lowPart);
    }
    static constexpr fixed_width_uint get_difference(fixed_width_uint a, fixed_width_uint b)
    {
        return fixed_width_uint(get_difference_borrow(a.lowPart, b.lowPart) ? a.highPart - b.highPart - (half_width_type)1 : a.highPart - b.highPart, a.lowPart - b.lowPart);
    }
#if 1 // use new division algorithm
    static cpp14_constexpr fixed_width_uint divide(fixed_width_uint dividend, fixed_width_uint divisor)
    {
        return div_imp(dividend, divisor);
    }
    static cpp14_constexpr std::pair<fixed_width_uint, fixed_width_uint> divmod(fixed_width_uint dividend, fixed_width_uint divisor)
    {
        fixed_width_uint quotient = divide(dividend, divisor);
        return std::pair<fixed_width_uint, fixed_width_uint>(quotient, dividend - quotient * divisor);
    }
    static cpp14_constexpr fixed_width_uint modulus(fixed_width_uint dividend, fixed_width_uint divisor)
    {
        return std::get<1>(divmod(dividend, divisor));
    }
#else
    static cpp14_constexpr std::pair<fixed_width_uint, fixed_width_uint> divmod_helper(fixed_width_uint dividend, fixed_width_uint divisor, int bit_index, fixed_width_uint quotient)
    {
        for(; bit_index >= 0; bit_index--)
        {
            if(dividend >= (divisor << bit_index))
            {
                dividend -= (divisor << bit_index);
                quotient |= (fixed_width_uint(1) << bit_index);
            }
        }
        return std::pair<fixed_width_uint, fixed_width_uint>(quotient, dividend);
    }
    static cpp14_constexpr std::pair<fixed_width_uint, fixed_width_uint> divmod(fixed_width_uint dividend, fixed_width_uint divisor)
    {
        return dividend.highPart == (half_width_type)0 && divisor.highPart == (half_width_type)0 ? std::pair<fixed_width_uint, fixed_width_uint>(fixed_width_uint(half_width_type(0), dividend.lowPart / divisor.lowPart), fixed_width_uint(half_width_type(0), dividend.lowPart % divisor.lowPart)) :
            dividend == divisor ? std::pair<fixed_width_uint, fixed_width_uint>(fixed_width_uint(1), fixed_width_uint(0)) :
            dividend < divisor ? std::pair<fixed_width_uint, fixed_width_uint>(fixed_width_uint(0), divisor) :
            divmod_helper(dividend, divisor, ilog_base_2(dividend) - ilog_base_2(divisor), fixed_width_uint(0));
    }
    static cpp14_constexpr fixed_width_uint divide(fixed_width_uint dividend, fixed_width_uint divisor)
    {
        return std::get<0>(divmod(dividend, divisor));
    }
    static cpp14_constexpr fixed_width_uint modulus(fixed_width_uint dividend, fixed_width_uint divisor)
    {
        return std::get<1>(divmod(dividend, divisor));
    }
#endif
    template <size_t N, typename = typename std::enable_if<(N > bit_count) && (N & (N - 1)) == 0>::type>
    static constexpr fixed_width_uint convert_size(fixed_width_uint<N> v)
    {
        return convert_size(v.lowPart);
    }
    static constexpr fixed_width_uint convert_size(fixed_width_uint v)
    {
        return v;
    }
    template <size_t N, typename = void, typename = typename std::enable_if<(N < bit_count) && (N & (N - 1)) == 0>::type>
    static constexpr fixed_width_uint convert_size(fixed_width_uint<N> v)
    {
        return convert_size(fixed_width_uint<N * 2>((fixed_width_uint<N>)0, v));
    }
public:
    fixed_width_uint()
    {
    }
    explicit constexpr fixed_width_uint(std::uintmax_t value)
        : lowPart(bit_count >= 2 * uintmax_bit_count ? value : value & (((std::uintmax_t)1 << bit_count / 2) - 1U)), highPart(bit_count >= 2 * uintmax_bit_count ? 0 : (value >> (bit_count / 2)) & (((std::uintmax_t)1 << bit_count / 2) - 1U))
    {
    }
    explicit constexpr operator std::uintmax_t() const
    {
        return bit_count >= 2 * uintmax_bit_count ? (std::uintmax_t)lowPart : ((std::uintmax_t)highPart << half_width_type::bit_count) | (std::uintmax_t)lowPart;
    }
    template <size_t N>
    explicit constexpr fixed_width_uint(fixed_width_uint<N> rt)
        : fixed_width_uint(convert_size(rt))
    {
    }
    friend constexpr fixed_width_uint operator +(fixed_width_uint a, fixed_width_uint b)
    {
        return get_sum(a, b);
    }
    friend constexpr fixed_width_uint operator -(fixed_width_uint a, fixed_width_uint b)
    {
        return get_difference(a, b);
    }
    constexpr fixed_width_uint operator +() const
    {
        return *this;
    }
    constexpr fixed_width_uint operator -() const
    {
        return lowPart == (half_width_type)0 ? fixed_width_uint(-highPart, (half_width_type)0) : fixed_width_uint(~highPart, -lowPart);
    }
    constexpr fixed_width_uint operator ~() const
    {
        return fixed_width_uint(~highPart, ~lowPart);
    }
    cpp14_constexpr const fixed_width_uint &operator ++()
    {
        if(++lowPart == (half_width_type)0)
            ++highPart;
        return *this;
    }
    cpp14_constexpr const fixed_width_uint &operator --()
    {
        if(lowPart == (half_width_type)0)
            --highPart;
        --lowPart;
        return *this;
    }
    cpp14_constexpr fixed_width_uint operator ++(int)
    {
        fixed_width_uint retval = *this;
        operator ++();
        return retval;
    }
    cpp14_constexpr fixed_width_uint operator --(int)
    {
        fixed_width_uint retval = *this;
        operator --();
        return retval;
    }
    friend constexpr fixed_width_uint operator *(fixed_width_uint a, fixed_width_uint b)
    {
        return get_product(a, b);
    }
    friend constexpr fixed_width_uint operator &(fixed_width_uint a, fixed_width_uint b)
    {
        return fixed_width_uint(a.highPart & b.highPart, a.lowPart & b.lowPart);
    }
    friend constexpr fixed_width_uint operator |(fixed_width_uint a, fixed_width_uint b)
    {
        return fixed_width_uint(a.highPart | b.highPart, a.lowPart | b.lowPart);
    }
    friend constexpr fixed_width_uint operator ^(fixed_width_uint a, fixed_width_uint b)
    {
        return fixed_width_uint(a.highPart ^ b.highPart, a.lowPart ^ b.lowPart);
    }
    cpp14_constexpr const fixed_width_uint &operator +=(fixed_width_uint r)
    {
        return *this = get_sum(*this, r);
    }
    cpp14_constexpr const fixed_width_uint &operator -=(fixed_width_uint r)
    {
        return *this = get_difference(*this, r);
    }
    cpp14_constexpr const fixed_width_uint &operator *=(fixed_width_uint r)
    {
        return *this = get_product(*this, r);
    }
    cpp14_constexpr const fixed_width_uint &operator &=(fixed_width_uint r)
    {
        highPart &= r.highPart;
        lowPart &= r.lowPart;
        return *this;
    }
    cpp14_constexpr const fixed_width_uint &operator |=(fixed_width_uint r)
    {
        highPart |= r.highPart;
        lowPart |= r.lowPart;
        return *this;
    }
    cpp14_constexpr const fixed_width_uint &operator ^=(fixed_width_uint r)
    {
        highPart ^= r.highPart;
        lowPart ^= r.lowPart;
        return *this;
    }
    friend constexpr fixed_width_uint operator <<(fixed_width_uint a, size_t b)
    {
        return b >= bit_count ? fixed_width_uint(0) :
            b == 0 ? a :
            b == half_width_type::bit_count ? fixed_width_uint(a.lowPart, half_width_type(0)) :
            b > half_width_type::bit_count ? fixed_width_uint(a.lowPart << (b - half_width_type::bit_count), half_width_type(0)) :
            fixed_width_uint((a.highPart << b) | (a.lowPart >> (half_width_type::bit_count - b)), a.lowPart << b);
    }
    friend constexpr fixed_width_uint operator >>(fixed_width_uint a, size_t b)
    {
        return b >= bit_count ? fixed_width_uint(0) :
            b == 0 ? a :
            b == half_width_type::bit_count ? fixed_width_uint(half_width_type(0), a.highPart) :
            b > half_width_type::bit_count ? fixed_width_uint(half_width_type(0), a.highPart >> (b - half_width_type::bit_count)) :
            fixed_width_uint(a.highPart >> b, (a.highPart << (half_width_type::bit_count - b)) | (a.lowPart >> b));
    }
    cpp14_constexpr const fixed_width_uint &operator <<=(size_t r)
    {
        return *this = *this << r;
    }
    cpp14_constexpr const fixed_width_uint &operator >>=(size_t r)
    {
        return *this = *this >> r;
    }
    friend constexpr bool operator ==(fixed_width_uint a, fixed_width_uint b)
    {
        return a.highPart == b.highPart && a.lowPart == b.lowPart;
    }
    friend constexpr bool operator !=(fixed_width_uint a, fixed_width_uint b)
    {
        return !(a == b);
    }
    friend constexpr bool operator <(fixed_width_uint a, fixed_width_uint b)
    {
        return a.highPart < b.highPart || (a.highPart == b.highPart && a.lowPart < b.lowPart);
    }
    friend constexpr bool operator >=(fixed_width_uint a, fixed_width_uint b)
    {
        return !(a < b);
    }
    friend constexpr bool operator <=(fixed_width_uint a, fixed_width_uint b)
    {
        return !(a > b);
    }
    friend constexpr bool operator >(fixed_width_uint a, fixed_width_uint b)
    {
        return b < a;
    }
    friend cpp14_constexpr fixed_width_uint operator /(fixed_width_uint a, fixed_width_uint b)
    {
        return divide(a, b);
    }
    friend cpp14_constexpr fixed_width_uint operator %(fixed_width_uint a, fixed_width_uint b)
    {
        return modulus(a, b);
    }
    cpp14_constexpr const fixed_width_uint &operator /=(fixed_width_uint r)
    {
        return *this = *this / r;
    }
    cpp14_constexpr const fixed_width_uint &operator %=(fixed_width_uint r)
    {
        return *this = *this % r;
    }
    friend std::ostream &operator <<(std::ostream &os, fixed_width_uint value)
    {
        if(!os)
            return os;
        if(bit_count <= uintmax_bit_count)
        {
            return os << (uintmax_t)value;
        }
        std::string str = "";
        fixed_width_uint chunk_modulus = fixed_width_uint((uintmax_t)1000000000 * (uintmax_t)100000 * (uintmax_t)100000);
        size_t chunk_digit_count = 19;
        bool can_use_shift = false;
        size_t chunk_shift_count = 0;
        const char *base_prefix = "";
        switch(os.flags() & std::ios_base::basefield)
        {
        case std::ios_base::oct:
            chunk_digit_count = 21;
            chunk_shift_count = 3 * chunk_digit_count;
            can_use_shift = true;
            chunk_modulus = fixed_width_uint(1) << chunk_shift_count;
            base_prefix = "0";
            break;
        case std::ios_base::hex:
            chunk_digit_count = 16;
            chunk_shift_count = 4 * chunk_digit_count;
            can_use_shift = true;
            chunk_modulus = fixed_width_uint(1) << chunk_shift_count;
            base_prefix = "0x";
            if(os.flags() & std::ios_base::uppercase)
                base_prefix = "0X";
            break;
        case std::ios_base::dec:
        default:
            break;
        }

        while(value >= chunk_modulus)
        {
            uintmax_t chunk_value;
            if(can_use_shift)
            {
                chunk_value = (uintmax_t)((chunk_modulus - (fixed_width_uint)1) & value);
                value >>= chunk_shift_count;
            }
            else
            {
                chunk_value = (uint64_t)(value % chunk_modulus);
                value /= chunk_modulus;
            }
            std::ostringstream ss;
            std::ios_base::fmtflags flags_mask = std::ios_base::uppercase | std::ios_base::basefield;
            ss.flags((flags_mask & os.flags()) | (ss.flags() & ~flags_mask));
            ss << chunk_value;
            std::string chunk_str = ss.str();
            if(chunk_str.size() < chunk_digit_count)
                chunk_str.insert(0, chunk_digit_count - chunk_str.size(), '0');
            str = chunk_str + str;
        }
        std::ostringstream ss;
        std::ios_base::fmtflags flags_mask = std::ios_base::uppercase | std::ios_base::basefield;
        ss.flags((flags_mask & os.flags()) | (ss.flags() & ~flags_mask));
        ss << (uintmax_t)value;
        str = ss.str() + str;

        if(os.flags() & std::ios_base::showbase)
            str = base_prefix + str;

        if(os.flags() & std::ios_base::showpos && (os.flags() & std::ios_base::basefield) == std::ios_base::dec)
            str = "+" + str;

        size_t width = os.width() < 0 ? (size_t)0 : (size_t)os.width();

        if(width > str.size())
        {
            if((os.flags() & std::ios_base::adjustfield) == std::ios_base::left)
                str.append(width - str.size(), os.fill());
            else
                str.insert(0, width - str.size(), os.fill());
        }

        os << str;
        os.width(0);
        return os;
    }
private:
    static constexpr fixed_width_uint floor_sqrt_helper(fixed_width_uint v, fixed_width_uint bit, fixed_width_uint retval)
    {
        // use new sqrt algorithm once we speed up division
#if __cplusplus > 201103L
        for(; bit != (fixed_width_uint)0; bit = bit >> 2)
        {
            if(v >= bit + retval)
            {
                v = v - (bit + retval);
                retval = (retval >> 1) | bit;
            }
            else
            {
                retval = retval >> 1;
            }
        }
        return retval;
#else
        return bit == (fixed_width_uint)0 ? retval :
            v >= bit + retval ? floor_sqrt_helper(v - (bit + retval), bit >> 2, (retval >> 1) | bit) :
            floor_sqrt_helper(v, bit >> 2, retval >> 1);
#endif
    }
public:
    friend constexpr fixed_width_uint floor_sqrt(fixed_width_uint v)
    {
        return floor_sqrt_helper(v, (fixed_width_uint)1 << ((bit_count - 1) & ~1), (fixed_width_uint)0);
    }
    friend constexpr intmax_t ilog2(fixed_width_uint v)
    {
        return v.highPart == (half_width_type)0 ? ilog2(v.lowPart) : half_width_type::bit_count + ilog2(v.highPart);
    }
    friend constexpr fixed_width_uint<bit_count * 2> wide_multiply(fixed_width_uint a, fixed_width_uint b)
    {
        return fixed_width_uint_multiply_helper<bit_count>::multiply(a, b);
    }
    static cpp14_constexpr fixed_width_uint multiplicative_inverse(fixed_width_uint d_in);
private:
    static cpp14_constexpr fixed_width_uint div_imp(fixed_width_uint n_in, fixed_width_uint d_in);
};

template <size_t bitCount>
constexpr typename fixed_width_uint<bitCount>::half_width_type fixed_width_uint<bitCount>::quarter_mask;

template <size_t bitCount>
struct fixed_width_int final
{
    template <size_t>
    friend struct fixed_width_int;
    static constexpr size_t bit_count = bitCount;
    static_assert(bit_count >= 8, "invalid bit count for fixed_width_int : must be at least 8");
    static_assert((bit_count & (bit_count - 1)) == 0, "invalid bit count for fixed_width_int : must be a power of 2");
private:
    static constexpr fixed_width_uint<bit_count> sign_bit = fixed_width_uint<bit_count>(1) << (bit_count - 1);
    fixed_width_uint<bit_count> value;
public:
    template <size_t N>
    explicit constexpr fixed_width_int(fixed_width_uint<N> value)
        : value(value)
    {
    }
    fixed_width_int()
    {
    }
    explicit constexpr fixed_width_int(intmax_t value)
        : value(value < 0 ? -fixed_width_uint<bit_count>((uintmax_t)-value) : fixed_width_uint<bit_count>((uintmax_t)value))
    {
    }
    template <size_t N, typename = typename std::enable_if<(N > bit_count)>::type>
    explicit constexpr fixed_width_int(fixed_width_int<N> value)
        : value(value.value)
    {
    }
    template <size_t N, typename = void, typename = typename std::enable_if<(N < bit_count)>::type>
    explicit constexpr fixed_width_int(fixed_width_int<N> value)
        : value(value < (fixed_width_int<N>)0 ? -(fixed_width_uint<bit_count>)-value : (fixed_width_uint<bit_count>)value)
    {
    }
    template <size_t N>
    explicit constexpr operator fixed_width_uint<N>() const
    {
        return ((fixed_width_int<N>)*this).value;
    }
    explicit constexpr operator intmax_t() const
    {
        return (value & sign_bit) != (fixed_width_uint<bit_count>)0 ? -(intmax_t)(uintmax_t)-value : (intmax_t)(uintmax_t)value;
    }
    friend constexpr fixed_width_int operator +(fixed_width_int a, fixed_width_int b)
    {
        return fixed_width_int(a.value + b.value);
    }
    friend constexpr fixed_width_int operator -(fixed_width_int a, fixed_width_int b)
    {
        return fixed_width_int(a.value - b.value);
    }
    friend constexpr fixed_width_int operator *(fixed_width_int a, fixed_width_int b)
    {
        return fixed_width_int(a.value * b.value);
    }
    friend constexpr fixed_width_int operator &(fixed_width_int a, fixed_width_int b)
    {
        return fixed_width_int(a.value & b.value);
    }
    friend constexpr fixed_width_int operator |(fixed_width_int a, fixed_width_int b)
    {
        return fixed_width_int(a.value | b.value);
    }
    friend constexpr fixed_width_int operator ^(fixed_width_int a, fixed_width_int b)
    {
        return fixed_width_int(a.value ^ b.value);
    }
    friend constexpr fixed_width_int operator <<(fixed_width_int a, size_t b)
    {
        return fixed_width_int(a.value << b);
    }
    friend constexpr fixed_width_int operator >>(fixed_width_int a, size_t b)
    {
        return fixed_width_int((a.value & sign_bit) != (fixed_width_uint<bit_count>)0 ? ~(~a.value >> b) : a.value >> b);
    }
    friend constexpr fixed_width_int operator /(fixed_width_int a, fixed_width_int b)
    {
        return (a.value & sign_bit) != (fixed_width_uint<bit_count>)0 ? ((b.value & sign_bit) != (fixed_width_uint<bit_count>)0 ? fixed_width_int(-a.value / -b.value) : fixed_width_int(-(-a.value / b.value))) : ((b.value & sign_bit) != (fixed_width_uint<bit_count>)0 ? fixed_width_int(-(a.value / -b.value)) : fixed_width_int(a.value / b.value));
    }
    friend constexpr fixed_width_int operator %(fixed_width_int a, fixed_width_int b)
    {
        return (a.value & sign_bit) != (fixed_width_uint<bit_count>)0 ? ((b.value & sign_bit) != (fixed_width_uint<bit_count>)0 ? fixed_width_int(-(-a.value % -b.value)) : fixed_width_int(-(-a.value % b.value))) : ((b.value & sign_bit) != (fixed_width_uint<bit_count>)0 ? fixed_width_int(a.value / -b.value) : fixed_width_int(a.value / b.value));
    }
    cpp14_constexpr const fixed_width_int &operator ++()
    {
        ++value;
        return *this;
    }
    cpp14_constexpr const fixed_width_int &operator --()
    {
        --value;
        return *this;
    }
    cpp14_constexpr fixed_width_int operator ++(int)
    {
        return fixed_width_int(value++);
    }
    cpp14_constexpr fixed_width_int operator --(int)
    {
        return fixed_width_int(value--);
    }
    constexpr fixed_width_int operator -() const
    {
        return fixed_width_int(-value);
    }
    constexpr fixed_width_int operator ~() const
    {
        return fixed_width_int(~value);
    }
    constexpr const fixed_width_int &operator +() const
    {
        return *this;
    }
    cpp14_constexpr const fixed_width_int &operator +=(fixed_width_int r)
    {
        return *this = *this + r;
    }
    cpp14_constexpr const fixed_width_int &operator -=(fixed_width_int r)
    {
        return *this = *this - r;
    }
    cpp14_constexpr const fixed_width_int &operator *=(fixed_width_int r)
    {
        return *this = *this * r;
    }
    cpp14_constexpr const fixed_width_int &operator /=(fixed_width_int r)
    {
        return *this = *this / r;
    }
    cpp14_constexpr const fixed_width_int &operator %=(fixed_width_int r)
    {
        return *this = *this % r;
    }
    cpp14_constexpr const fixed_width_int &operator &=(fixed_width_int r)
    {
        return *this = *this & r;
    }
    cpp14_constexpr const fixed_width_int &operator |=(fixed_width_int r)
    {
        return *this = *this | r;
    }
    cpp14_constexpr const fixed_width_int &operator ^=(fixed_width_int r)
    {
        return *this = *this ^ r;
    }
    cpp14_constexpr const fixed_width_int &operator <<=(size_t r)
    {
        return *this = *this << r;
    }
    cpp14_constexpr const fixed_width_int &operator >>=(size_t r)
    {
        return *this = *this >> r;
    }
    friend constexpr bool operator ==(fixed_width_int a, fixed_width_int b)
    {
        return a.value == b.value;
    }
    friend constexpr bool operator !=(fixed_width_int a, fixed_width_int b)
    {
        return a.value != b.value;
    }
    friend constexpr bool operator <(fixed_width_int a, fixed_width_int b)
    {
        return (a.value ^ sign_bit) < (b.value ^ sign_bit);
    }
    friend constexpr bool operator >(fixed_width_int a, fixed_width_int b)
    {
        return (a.value ^ sign_bit) > (b.value ^ sign_bit);
    }
    friend constexpr bool operator <=(fixed_width_int a, fixed_width_int b)
    {
        return (a.value ^ sign_bit) <= (b.value ^ sign_bit);
    }
    friend constexpr bool operator >=(fixed_width_int a, fixed_width_int b)
    {
        return (a.value ^ sign_bit) >= (b.value ^ sign_bit);
    }
    friend std::ostream &operator <<(std::ostream &os, fixed_width_int value)
    {
        if(!os)
            return os;
        std::ostringstream ss;
        std::ios_base::fmtflags flags_mask = std::ios_base::uppercase | std::ios_base::basefield;
        ss.flags((flags_mask & os.flags()) | (ss.flags() & ~flags_mask));
        bool is_neg = (value.value & sign_bit) == sign_bit, can_show_sign = true;
        if((os.flags() & std::ios_base::basefield) != std::ios_base::dec)
        {
            is_neg = false;
            can_show_sign = false;
        }
        if(is_neg)
            value = -value;
        ss << value.value;
        std::string str = ss.str();
        const char *base_prefix = "";
        switch(os.flags() & std::ios_base::basefield)
        {
        case std::ios_base::oct:
            base_prefix = "0";
            break;
        case std::ios_base::hex:
            base_prefix = "0x";
            if(os.flags() & std::ios_base::uppercase)
                base_prefix = "0X";
            break;
        case std::ios_base::dec:
        default:
            break;
        }

        if(os.flags() & std::ios_base::showbase)
            str = base_prefix + str;

        if(is_neg)
            str = "-" + str;
        else if(os.flags() & std::ios_base::showpos && can_show_sign)
            str = "+" + str;

        size_t width = os.width() < 0 ? (size_t)0 : (size_t)os.width();

        if(width > str.size())
        {
            if((os.flags() & std::ios_base::adjustfield) == std::ios_base::left)
                str.append(width - str.size(), os.fill());
            else
                str.insert(0, width - str.size(), os.fill());
        }

        os << str;
        os.width(0);
        return os;
    }
    static constexpr fixed_width_int abs(fixed_width_int v)
    {
        return (v.value & sign_bit) == (fixed_width_uint<bit_count>)0 ? v : -v;
    }
    friend constexpr intmax_t ilog2(fixed_width_int v)
    {
        return v <= (fixed_width_int)0 ? -1 : ilog2((fixed_width_uint<bit_count>)v);
    }
    friend constexpr fixed_width_int<bit_count * 2> wide_multiply(fixed_width_int a, fixed_width_int b)
    {
        return a >= (fixed_width_int)0 ?
            (b >= (fixed_width_int)0 ?
                (fixed_width_int<bit_count * 2>)wide_multiply((fixed_width_uint<bit_count>)a, (fixed_width_uint<bit_count>)b) :
                -(fixed_width_int<bit_count * 2>)wide_multiply((fixed_width_uint<bit_count>)a, (fixed_width_uint<bit_count>)-b)) :
            (b >= (fixed_width_int)0 ?
                -(fixed_width_int<bit_count * 2>)wide_multiply((fixed_width_uint<bit_count>)-a, (fixed_width_uint<bit_count>)b) :
                (fixed_width_int<bit_count * 2>)wide_multiply((fixed_width_uint<bit_count>)-a, (fixed_width_uint<bit_count>)-b));
    }
};

template <size_t N>
constexpr fixed_width_int<N> abs(fixed_width_int<N> v)
{
    return fixed_width_int<N>::abs(v);
}

template <size_t bitCount>
constexpr fixed_width_uint<fixed_width_int<bitCount>::bit_count> fixed_width_int<bitCount>::sign_bit;

namespace std
{
template <size_t N>
struct numeric_limits<fixed_width_uint<N>>
{
private:
    static constexpr int constexpr_ifloor(long double v)
    {
        return v < (int)v ? (int)v - 1 : (int)v;
    }
    static constexpr int constexpr_iceil(long double v)
    {
        return v > (int)v ? (int)v + 1 : (int)v;
    }
    static constexpr long double log10_2 = 0.30102999566398119521373889472449302676818988L;
public:
    static constexpr bool is_specialized = true;
    static constexpr fixed_width_uint<N> min()
    {
        return fixed_width_uint<N>(0);
    }
    static constexpr fixed_width_uint<N> max()
    {
        return ~fixed_width_uint<N>(0);
    }
    static constexpr fixed_width_uint<N> lowest()
    {
        return fixed_width_uint<N>(0);
    }
    static constexpr int digits = fixed_width_uint<N>::bit_count;
    static constexpr int digits10 = constexpr_ifloor(log10_2 * digits);
    static constexpr int max_digits10 = constexpr_iceil(log10_2 * digits);
    static constexpr bool is_signed = false;
    static constexpr bool is_integer = true;
    static constexpr bool is_exact = true;
    static constexpr int radix = 2;
    static constexpr fixed_width_uint<N> epsilon()
    {
        return fixed_width_uint<N>(1);
    }
    static constexpr fixed_width_uint<N> round_error()
    {
        return fixed_width_uint<N>(0);
    }
    static constexpr int min_exponent = 0;
    static constexpr int min_exponent10 = 0;
    static constexpr int max_exponent = digits;
    static constexpr int max_exponent10 = digits10;
    static constexpr bool has_infinity = false;
    static constexpr bool has_quiet_NaN = false;
    static constexpr bool has_signaling_NaN = false;
    static constexpr float_denorm_style has_denorm = float_denorm_style::denorm_absent;
    static constexpr bool has_denorm_loss = false;
    static constexpr fixed_width_uint<N> infinity()
    {
        return max();
    }
    static constexpr fixed_width_uint<N> quiet_NaN()
    {
        return fixed_width_uint<N>(0);
    }
    static constexpr fixed_width_uint<N> signaling_NaN()
    {
        return fixed_width_uint<N>(0);
    }
    static constexpr fixed_width_uint<N> denorm_min()
    {
        return min();
    }
    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = true;
    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;
    static constexpr float_round_style round_style = float_round_style::round_toward_zero;
};
template <size_t N>
struct numeric_limits<const fixed_width_uint<N>> : public numeric_limits<fixed_width_uint<N>>
{
};

template <size_t N>
struct numeric_limits<fixed_width_int<N>>
{
private:
    static constexpr int constexpr_ifloor(long double v)
    {
        return v < (int)v ? (int)v - 1 : (int)v;
    }
    static constexpr int constexpr_iceil(long double v)
    {
        return v > (int)v ? (int)v + 1 : (int)v;
    }
    static constexpr long double log10_2 = 0.30102999566398119521373889472449302676818988L;
public:
    static constexpr bool is_specialized = true;
    static constexpr fixed_width_int<N> min()
    {
        return fixed_width_int<N>(1) << (fixed_width_int<N>::bit_count - 1);
    }
    static constexpr fixed_width_int<N> max()
    {
        return ~(fixed_width_int<N>(1) << (fixed_width_int<N>::bit_count - 1));
    }
    static constexpr fixed_width_int<N> lowest()
    {
        return min();
    }
    static constexpr int digits = fixed_width_int<N>::bit_count - 1;
    static constexpr int digits10 = constexpr_ifloor(log10_2 * digits);
    static constexpr int max_digits10 = constexpr_iceil(log10_2 * digits);
    static constexpr bool is_signed = true;
    static constexpr bool is_integer = true;
    static constexpr bool is_exact = true;
    static constexpr int radix = 2;
    static constexpr fixed_width_int<N> epsilon()
    {
        return fixed_width_int<N>(1);
    }
    static constexpr fixed_width_int<N> round_error()
    {
        return fixed_width_int<N>(0);
    }
    static constexpr int min_exponent = 0;
    static constexpr int min_exponent10 = 0;
    static constexpr int max_exponent = digits;
    static constexpr int max_exponent10 = digits10;
    static constexpr bool has_infinity = false;
    static constexpr bool has_quiet_NaN = false;
    static constexpr bool has_signaling_NaN = false;
    static constexpr float_denorm_style has_denorm = float_denorm_style::denorm_absent;
    static constexpr bool has_denorm_loss = false;
    static constexpr fixed_width_int<N> infinity()
    {
        return max();
    }
    static constexpr fixed_width_int<N> quiet_NaN()
    {
        return fixed_width_int<N>(0);
    }
    static constexpr fixed_width_int<N> signaling_NaN()
    {
        return fixed_width_int<N>(0);
    }
    static constexpr fixed_width_int<N> denorm_min()
    {
        return min();
    }
    static constexpr bool is_iec559 = false;
    static constexpr bool is_bounded = true;
    static constexpr bool is_modulo = true;
    static constexpr bool traps = false;
    static constexpr bool tinyness_before = false;
    static constexpr float_round_style round_style = float_round_style::round_toward_zero;
};
template <size_t N>
struct numeric_limits<const fixed_width_int<N>> : public numeric_limits<fixed_width_int<N>>
{
};
}

template <typename T, size_t bitCount, typename childClass, typename biggestIntType>
class builtin_fixed_width_int
{
    template <size_t N>
    friend struct fixed_width_int;
    template <size_t N>
    friend struct fixed_width_uint;
    T value;
protected:
    struct make_flag_t
    {
    };
    constexpr builtin_fixed_width_int(T value, make_flag_t)
        : value(value)
    {
    }
    static constexpr childClass make(T value)
    {
        return childClass(value, make_flag_t());
    }
public:
    static constexpr size_t bit_count = bitCount;
    builtin_fixed_width_int()
    {
    }
    explicit constexpr builtin_fixed_width_int(biggestIntType v)
        : value(v)
    {
    }
    explicit constexpr operator biggestIntType() const
    {
        return value;
    }
    constexpr childClass operator +() const
    {
        return make(value);
    }
    constexpr childClass operator -() const
    {
        return make(-value);
    }
    constexpr childClass operator ~() const
    {
        return make(~value);
    }
    cpp14_constexpr childClass operator ++()
    {
        return make(++value);
    }
    cpp14_constexpr childClass operator --()
    {
        return make(--value);
    }
    cpp14_constexpr childClass operator ++(int)
    {
        return make(value++);
    }
    cpp14_constexpr childClass operator --(int)
    {
        return make(value--);
    }
    friend constexpr childClass operator *(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return make(a.value * b.value);
    }
    friend constexpr childClass operator /(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return make(a.value / b.value);
    }
    friend constexpr childClass operator %(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return make(a.value % b.value);
    }
    friend constexpr childClass operator +(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return make(a.value + b.value);
    }
    friend constexpr childClass operator -(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return make(a.value - b.value);
    }
    friend constexpr childClass operator <<(builtin_fixed_width_int a, int b)
    {
        return make(a.value << b);
    }
    friend constexpr childClass operator >>(builtin_fixed_width_int a, int b)
    {
        return make(a.value >> b);
    }
    friend constexpr bool operator <(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return a.value < b.value;
    }
    friend constexpr bool operator >(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return a.value > b.value;
    }
    friend constexpr bool operator <=(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return a.value <= b.value;
    }
    friend constexpr bool operator >=(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return a.value >= b.value;
    }
    friend constexpr bool operator ==(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return a.value == b.value;
    }
    friend constexpr bool operator !=(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return a.value != b.value;
    }
    friend constexpr childClass operator &(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return make(a.value & b.value);
    }
    friend constexpr childClass operator |(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return make(a.value | b.value);
    }
    friend constexpr childClass operator ^(builtin_fixed_width_int a, builtin_fixed_width_int b)
    {
        return make(a.value ^ b.value);
    }
    cpp14_constexpr childClass operator *=(builtin_fixed_width_int r)
    {
        return make(value *= r.value);
    }
    cpp14_constexpr childClass operator /=(builtin_fixed_width_int r)
    {
        return make(value /= r.value);
    }
    cpp14_constexpr childClass operator %=(builtin_fixed_width_int r)
    {
        return make(value %= r.value);
    }
    cpp14_constexpr childClass operator +=(builtin_fixed_width_int r)
    {
        return make(value += r.value);
    }
    cpp14_constexpr childClass operator -=(builtin_fixed_width_int r)
    {
        return make(value -= r.value);
    }
    cpp14_constexpr childClass operator <<=(size_t r)
    {
        return make(value <<= r);
    }
    cpp14_constexpr childClass operator >>=(size_t r)
    {
        return make(value >>= r);
    }
    cpp14_constexpr childClass operator &=(builtin_fixed_width_int r)
    {
        return make(value &= r.value);
    }
    cpp14_constexpr childClass operator |=(builtin_fixed_width_int r)
    {
        return make(value |= r.value);
    }
    cpp14_constexpr childClass operator ^=(builtin_fixed_width_int r)
    {
        return make(value ^= r.value);
    }
    static constexpr childClass abs(builtin_fixed_width_int v)
    {
        return make(v.value < 0 ? -v.value : v.value);
    }
    friend std::ostream &operator <<(std::ostream &os, builtin_fixed_width_int v)
    {
        return os << v.value;
    }
    friend std::istream &operator <<(std::istream &is, builtin_fixed_width_int &v)
    {
        return is >> v.value;
    }
private:
    static constexpr childClass floor_sqrt_helper(builtin_fixed_width_int v, builtin_fixed_width_int bit, childClass retval)
    {
        return bit.value == 0 ? retval :
            v.value >= bit.value + retval.value ? floor_sqrt_helper(v - (bit + retval), bit >> 2, (retval >> 1) | bit) :
            floor_sqrt_helper(v, bit >> 2, retval >> 1);
    }
    static constexpr intmax_t ilog2(uintmax_t v)
    {
        return v <= 0 ? -1 :
            (v >> 32) != 0 ? 32 + ilog2(v >> 32) :
            (v >> 16) != 0 ? 16 + ilog2(v >> 16) :
            (v >> 8) != 0 ? 8 + ilog2(v >> 8) :
            (v >> 4) != 0 ? 4 + ilog2(v >> 4) :
            (v >> 2) != 0 ? 2 + ilog2(v >> 2) :
            (v >> 1) != 0 ? 1 + ilog2(v >> 1) :
            0;
    }
public:
    friend constexpr childClass floor_sqrt(builtin_fixed_width_int v)
    {
        return floor_sqrt_helper(v, make(1) << ((bit_count - 1) & ~1), make(0));
    }
    friend constexpr intmax_t ilog2(builtin_fixed_width_int v)
    {
        return v.value <= 0 ? -1 :
            ilog2((uintmax_t)v.value);
    }
};

#define INSTANTIATE_FIXED_WIDTH_INT(n) \
template <> \
struct fixed_width_int<n> final : public builtin_fixed_width_int<std::int ## n ## _t, n, fixed_width_int<n>, std::intmax_t> \
{ \
    constexpr fixed_width_int(std::int ## n ## _t value, make_flag_t) \
        : builtin_fixed_width_int(value) \
    { \
    } \
    fixed_width_int() \
    { \
    } \
    explicit constexpr fixed_width_int(std::intmax_t v) \
        : builtin_fixed_width_int(v) \
    { \
    } \
    template <size_t N2> \
    explicit constexpr fixed_width_int(fixed_width_uint<N2> v) \
        : builtin_fixed_width_int((std::intmax_t)(std::uintmax_t)v)\
    { \
    } \
    template <size_t N2> \
    explicit constexpr fixed_width_int(fixed_width_int<N2> v) \
        : builtin_fixed_width_int((std::intmax_t)v)\
    { \
    } \
    template <typename T, typename = typename std::enable_if<std::is_same<fixed_width_int, T>::value>::type> \
    friend constexpr fixed_width_int<T::bit_count * 2> wide_multiply(fixed_width_int a, T b) \
    { \
        return fixed_width_int_multiply_helper<T::bit_count>::multiply(a, b); \
    } \
};

#define INSTANTIATE_FIXED_WIDTH_UINT(n) \
template <> \
struct fixed_width_uint<n> final : public builtin_fixed_width_int<std::uint ## n ## _t, n, fixed_width_uint<n>, std::uintmax_t> \
{ \
    constexpr fixed_width_uint(std::uint ## n ## _t value, make_flag_t) \
        : builtin_fixed_width_int(value) \
    { \
    } \
    fixed_width_uint() \
    { \
    } \
    explicit constexpr fixed_width_uint(std::uintmax_t v) \
        : builtin_fixed_width_int(v) \
    { \
    } \
    template <size_t N2> \
    explicit constexpr fixed_width_uint(fixed_width_uint<N2> v) \
        : builtin_fixed_width_int((std::uintmax_t)v)\
    { \
    } \
    template <size_t N2> \
    explicit constexpr fixed_width_uint(fixed_width_int<N2> v) \
        : builtin_fixed_width_int((std::uintmax_t)(std::intmax_t)v)\
    { \
    } \
    template <typename T, typename = typename std::enable_if<std::is_same<fixed_width_uint, T>::value>::type> \
    friend constexpr fixed_width_uint<T::bit_count * 2> wide_multiply(fixed_width_uint a, T b) \
    { \
        return fixed_width_uint_multiply_helper<T::bit_count>::multiply(a, b); \
    } \
};

// must instantiate from in ascending size order
// must instantiate uint's first
INSTANTIATE_FIXED_WIDTH_UINT(8)
INSTANTIATE_FIXED_WIDTH_INT(8)
INSTANTIATE_FIXED_WIDTH_UINT(16)
INSTANTIATE_FIXED_WIDTH_INT(16)
INSTANTIATE_FIXED_WIDTH_UINT(32)
INSTANTIATE_FIXED_WIDTH_INT(32)
INSTANTIATE_FIXED_WIDTH_UINT(64)
INSTANTIATE_FIXED_WIDTH_INT(64)

#undef INSTANTIATE_FIXED_WIDTH_INT
#undef INSTANTIATE_FIXED_WIDTH_UINT

#if 0
template <size_t bitCount>
struct fixed_width_uint_multiply_helper<bitCount, typename std::enable_if<(bitCount >= 64)>::type> // multiply using karatsuba's algorithm
{
    static constexpr size_t bit_count = bitCount;
private:
    static_assert(bit_count >= 16 && (bit_count & (bit_count - 1)) == 0, "invalid bit count");
    typedef fixed_width_uint<bit_count * 2> duint;
    typedef fixed_width_uint<bit_count> uint;
    typedef fixed_width_uint<bit_count / 2> huint;
    static constexpr std::pair<uint, uint> split_word_helper(uint v, uint high_part)
    {
        return std::pair<uint, uint>(high_part, v - (high_part << bit_count / 2));
    }
    static constexpr std::pair<uint, uint> split_word(uint v)
    {
        return split_word_helper(v, v >> bit_count / 2);
    }
    static constexpr uint multiply_part(uint a, uint b)
    {
        return a >= (uint)1 << huint::bit_count ?
            (b >= (uint)1 << huint::bit_count ?
                wide_multiply((huint)a, (huint)b) + ((a + b) << huint::bit_count) :
                wide_multiply((huint)a, (huint)b) + (b << huint::bit_count)) :
            (b >= (uint)1 << huint::bit_count ?
                wide_multiply((huint)a, (huint)b) + (a << huint::bit_count) :
                wide_multiply((huint)a, (huint)b));
    }
    static constexpr duint multiply_helper(uint ah, uint al, uint bh, uint bl, uint zl, uint zh)
    {
        return ((duint)zh << bit_count) + (duint)zl + ((duint)(multiply_part(al + ah, bl + bh) - zl - zh) << huint::bit_count);
    }
    static constexpr duint multiply_helper(uint ah, uint al, uint bh, uint bl)
    {
        return ah == (uint)0 && bh == (uint)0 ? (duint)multiply_part(al, bl) :
            multiply_helper(ah, al, bh, bl, multiply_part(al, bl), multiply_part(ah, bh));
    }
    static constexpr duint multiply_helper(std::pair<uint, uint> a, std::pair<uint, uint> b)
    {
        return multiply_helper(std::get<0>(a), std::get<1>(a), std::get<0>(b), std::get<1>(b));
    }
public:
    static constexpr duint multiply(uint a, uint b)
    {
        return multiply_helper(split_word(a), split_word(b));
    }
};
#else
template <size_t bitCount>
struct fixed_width_uint_multiply_helper<bitCount, typename std::enable_if<(bitCount >= 64)>::type>
{
    static constexpr size_t bit_count = bitCount;
private:
    static_assert(bit_count >= 16 && (bit_count & (bit_count - 1)) == 0, "invalid bit count");
    typedef fixed_width_uint<bit_count * 2> duint;
    typedef fixed_width_uint<bit_count> uint;
    typedef fixed_width_uint<bit_count / 2> huint;
    static constexpr std::pair<uint, uint> split_word_helper(uint v, uint high_part)
    {
        return std::pair<uint, uint>(high_part, v - (high_part << bit_count / 2));
    }
    static constexpr std::pair<uint, uint> split_word(uint v)
    {
        return split_word_helper(v, v >> bit_count / 2);
    }
    static constexpr uint multiply_part(uint a, uint b)
    {
        return wide_multiply((huint)a, (huint)b);
    }
    static constexpr duint multiply_helper(uint ah, uint al, uint bh, uint bl)
    {
        return ah == (uint)0 && bh == (uint)0 ? (duint)multiply_part(al, bl) :
            ((duint)multiply_part(ah, bh) << uint::bit_count) +
                (((duint)multiply_part(ah, bl) + (duint)multiply_part(al, bh)) << huint::bit_count) +
                (duint)multiply_part(al, bl);
    }
    static constexpr duint multiply_helper(std::pair<uint, uint> a, std::pair<uint, uint> b)
    {
        return multiply_helper(std::get<0>(a), std::get<1>(a), std::get<0>(b), std::get<1>(b));
    }
public:
    static constexpr duint multiply(uint a, uint b)
    {
        return multiply_helper(split_word(a), split_word(b));
    }
};
#endif

template <size_t bitCount>
struct fixed_width_uint_multiply_helper<bitCount, typename std::enable_if<(bitCount < 64)>::type> // multiply using built-in types
{
    static constexpr size_t bit_count = bitCount;
private:
    static_assert(bit_count >= 8 && (bit_count & (bit_count - 1)) == 0, "invalid bit count");
    typedef fixed_width_uint<bit_count * 2> duint;
    typedef fixed_width_uint<bit_count> uint;
public:
    static constexpr duint multiply(uint a, uint b)
    {
        return (duint)a * (duint)b;
    }
};

template <size_t bitCount>
struct fixed_width_int_multiply_helper<bitCount, typename std::enable_if<(bitCount >= 64)>::type> // multiply using karatsuba's algorithm
{
    static constexpr size_t bit_count = bitCount;
    static_assert(bit_count >= 16 && (bit_count & (bit_count - 1)) == 0, "invalid bit count");
    static constexpr fixed_width_int<bit_count * 2> multiply(fixed_width_int<bit_count> a, fixed_width_int<bit_count> b)
    {
        return a >= (fixed_width_int<bit_count>)0 ?
            (b >= (fixed_width_int<bit_count>)0 ?
                (fixed_width_int<bit_count * 2>)wide_multiply((fixed_width_uint<bit_count>)a, (fixed_width_uint<bit_count>)b) :
                -(fixed_width_int<bit_count * 2>)wide_multiply((fixed_width_uint<bit_count>)a, (fixed_width_uint<bit_count>)-b)) :
            (b >= (fixed_width_int<bit_count>)0 ?
                -(fixed_width_int<bit_count * 2>)wide_multiply((fixed_width_uint<bit_count>)-a, (fixed_width_uint<bit_count>)b) :
                (fixed_width_int<bit_count * 2>)wide_multiply((fixed_width_uint<bit_count>)-a, (fixed_width_uint<bit_count>)-b));
    }
};

template <size_t bitCount>
struct fixed_width_int_multiply_helper<bitCount, typename std::enable_if<(bitCount < 64)>::type> // multiply using built-in types
{
    static constexpr size_t bit_count = bitCount;
    static_assert(bit_count >= 8 && (bit_count & (bit_count - 1)) == 0, "invalid bit count");
    static constexpr fixed_width_int<bit_count * 2> multiply(fixed_width_int<bit_count> a, fixed_width_int<bit_count> b)
    {
        return (fixed_width_int<bit_count * 2>)a * (fixed_width_int<bit_count * 2>)b;
    }
};

template <size_t N, typename = typename std::enable_if<N == 8>::type>
static constexpr fixed_width_uint<8> make_repeated_byte(int v)
{
    return (fixed_width_uint<8>)(uintmax_t)v;
}

template <size_t N, typename = typename std::enable_if<(N > 8)>::type>
static constexpr fixed_width_uint<N> make_repeated_byte(int v)
{
    return make_repeated_byte_helper(make_repeated_byte<N / 2>(v));
}

template <size_t fraction_bits, size_t N>
static constexpr fixed_width_int<N> fractional_multiply(fixed_width_int<N> a, fixed_width_int<N> b)
{
    return (fixed_width_int<N>)(wide_multiply(a, b) >> fraction_bits);
}

template <size_t fraction_bits, size_t N>
static constexpr fixed_width_uint<N> fractional_multiply(fixed_width_uint<N> a, fixed_width_uint<N> b)
{
    return (fixed_width_uint<N>)(wide_multiply(a, b) >> fraction_bits);
}

template <size_t N>
static constexpr fixed_width_uint<N * 2> make_repeated_byte_helper(fixed_width_uint<N> v)
{
    return ((fixed_width_uint<N * 2>)v << N) + (fixed_width_uint<N * 2>)v;
}

template <size_t N>
cpp14_constexpr fixed_width_uint<N> fixed_width_uint<N>::multiplicative_inverse(fixed_width_uint<N> d_in)
{
    typedef fixed_width_uint<N> uint;
    typedef fixed_width_uint<N * 2> duint;
    typedef fixed_width_int<N * 2> dsint;
    constexpr size_t extra_precision = 4;
    if(d_in <= (uint)1)
        return ~(uint)0;
    if(d_in >= (uint)1 << (N - 1))
        return (uint)1;
    intmax_t shift_amount = N - ilog2(d_in) - 1;
    duint d = (duint)d_in << (shift_amount + extra_precision);
    constexpr duint const_48_17 = ((duint)make_repeated_byte<N>(0xD2) + ((duint)2 << N)) << extra_precision;
    constexpr duint const_32_17 = ((duint)make_repeated_byte<N>(0xE1) + ((duint)1 << N)) << extra_precision;
    constexpr dsint const_1 = (dsint)1 << (N + extra_precision);
    duint x = const_48_17 - fractional_multiply<N + extra_precision>(const_32_17, d);
    constexpr intmax_t iteration_count = 1 + ilog2((fixed_width_uint<64>)((N + 1) / 4)); // should be ceil(log2((N + 1)/log2(17))) but this is almost always the same and always >= the correct value
    for(intmax_t i = 0; i < iteration_count; i++)
    {
        auto d_x = fractional_multiply<N + extra_precision>((dsint)d, (dsint)x);
        auto e = const_1 - d_x;
        auto delta = fractional_multiply<N + extra_precision>((dsint)x, e);
        x = (duint)((dsint)x + delta);
    }
    shift_amount = N - shift_amount;
    x += (duint)1 << (shift_amount - 1 + extra_precision);
    return (uint)(x >> (shift_amount + extra_precision));
}

template <size_t N>
cpp14_constexpr fixed_width_uint<N> fixed_width_uint<N>::div_imp(fixed_width_uint<N> n_in, fixed_width_uint<N> d_in)
{
    typedef fixed_width_uint<N> uint;
    typedef fixed_width_uint<N * 2> wuint;
    if(d_in == (uint)0)
        throw std::domain_error("division by zero");
    if(d_in == (uint)1)
        return n_in;
    if(d_in == (uint)2)
        return n_in >> 1;
    if(d_in > n_in)
        return (uint)0;
    if((wuint)d_in << 1 > (wuint)n_in)
        return (uint)1;
    uint inv = multiplicative_inverse(d_in);
    wuint q = (wuint)fractional_multiply<N>(inv, n_in);
    q++;
    if((wuint)n_in < q * (wuint)d_in)
        q--;
    if((wuint)n_in < q * (wuint)d_in)
        q--;
    return (uint)q;
}

#endif // INT_TYPES_H_INCLUDED
