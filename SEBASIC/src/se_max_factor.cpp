#include "../include/se_max_factor.h"

/* The trial divisor increment wheel.  Use it to skip over divisors that
   are composites of 2, 3, 5, 7, or 11.  The part from WHEEL_START up to
   WHEEL_END is reused periodically, while the "lead in" is used to test
   for those primes and to jump onto the wheel.
   The first 4 elements correspond to the incremental offsets of the first 5
   primes (2 3 5 7 11).  The 5(th) element is the difference between that
   last prime and the next largest integer that is not a multiple of
   those primes.  The remaining numbers define the wheel.
   For more information, see

拭除增量轮盘.用它来跳过2 3 5 7 11的合数的除数.当"lead in"用来测试这些质数并在wheel上
跳跃，从WHEEL_START道WHEEL_END的部分被周期性重复使用.
头4个元素对应头五个质数（2 3 5 7 11）的增量。第5个元素是最后一个质数和不是这些质数的
倍数挨着的最大整数的差.剩下的数字定义了这个wheel.
   http://www.utm.edu/research/primes/glossary/WheelFactorization.html */
#define WHEEL_SIZE 5 
static const unsigned int wheel_tab[] = {
    1, 2, 2, 4, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4, 2, 6, 4,
    6, 8, 4, 2, 4, 2, 4, 14, 4, 6, 2, 10, 2, 6, 6, 4, 2, 4, 6, 2, 10, 2,
    4, 2, 12, 10, 2, 4, 2, 4, 6, 2, 6, 4, 6, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8,
    4, 2, 4, 6, 8, 6, 10, 2, 4, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 6, 10,
    2, 10, 2, 4, 2, 4, 6, 8, 4, 2, 4, 12, 2, 6, 4, 2, 6, 4, 6, 12, 2, 4,
    2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 10, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2,
    10, 2, 4, 6, 6, 2, 6, 6, 4, 6, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 6, 4,
    8, 6, 4, 6, 2, 4, 6, 8, 6, 4, 2, 10, 2, 6, 4, 2, 4, 2, 10, 2, 10, 2,
    4, 2, 4, 8, 6, 4, 2, 4, 6, 6, 2, 6, 4, 8, 4, 6, 8, 4, 2, 4, 2, 4, 8,
    6, 4, 6, 6, 6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2, 10, 2, 6,
    4, 6, 2, 6, 4, 2, 4, 6, 6, 8, 4, 2, 6, 10, 8, 4, 2, 4, 2, 4, 8, 10, 6,
    2, 4, 8, 6, 6, 4, 2, 4, 6, 2, 6, 4, 6, 2, 10, 2, 10, 2, 4, 2, 4, 6, 2,
    6, 4, 2, 4, 6, 6, 2, 6, 6, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 8, 4,
    6, 2, 6, 6, 4, 2, 4, 6, 8, 4, 2, 4, 2, 10, 2, 10, 2, 4, 2, 4, 6, 2,
    10, 2, 4, 6, 8, 6, 4, 2, 6, 4, 6, 8, 4, 6, 2, 4, 8, 6, 4, 6, 2, 4, 6,
    2, 6, 6, 4, 6, 6, 2, 6, 6, 4, 2, 10, 2, 10, 2, 4, 2, 4, 6, 2, 6, 4, 2,
    10, 6, 2, 6, 4, 2, 6, 4, 6, 8, 4, 2, 4, 2, 12, 6, 4, 6, 2, 4, 6, 2,
    12, 4, 2, 4, 8, 6, 4, 2, 4, 2, 10, 2, 10, 6, 2, 4, 6, 2, 6, 4, 2, 4,
    6, 6, 2, 6, 4, 2, 10, 6, 8, 6, 4, 2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 6,
    6, 4, 6, 2, 6, 4, 2, 4, 2, 10, 12, 2, 4, 2, 10, 2, 6, 4, 2, 4, 6, 6,
    2, 10, 2, 6, 4, 14, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2, 4, 6, 2, 6, 6, 4, 2,
    4, 6, 2, 6, 4, 2, 4, 12, 2, 12 
};

#define WHEEL_START (wheel_tab + WHEEL_SIZE)
#define WHEEL_END (wheel_tab + (sizeof wheel_tab / sizeof wheel_tab[0]))


int se_factor(uint64_t n0, uint64_t *factors)
{
    uint64_t n = n0, d, q;
    int n_factors = 0;
    unsigned int const *w = wheel_tab;

    if (n < 1)
        return n_factors;

    /* The exit condition in the following loop is correct because
       any time it is tested one of these 3 conditions holds:
       (1) d divides n
       (2) n is prime
       (3) n is composite but has no factors less than d.
       If (1) or (2) obviously the right thing happens.
       If (3), then since n is composite it is >= d^2. */
    d = 2;
    do {
        q = n / d;
        while (n == q * d) {
            factors[n_factors++] = d;    //后加
            n = q;
            q = n / d;
        }
        d += *(w++);
        if (w == WHEEL_END)       //循环
            w = WHEEL_START;
    } while (d <= q);  //<sqrt(n)

    if (n != 1 || n0 == 1) {
        factors[n_factors++] = n;
    }

    return n_factors;
}

//将number分成nfct个近似相等的因子
//数组 fct  个数 nfct
void se_split_factors(uint64_t number, uint64_t* fct, size_t nfct)/** Splits number in nfct aproximatively equal factors */
{
    uint64_t factors[SE_MAX_N_FACTORS];
    int i, ifct, n = se_factor(number, factors);
//初始化 nfct 已经固定大小
    for(ifct = 0; ifct < (int)nfct; ++ifct) {
        fct[ifct] = 1;
    }
//n 是 factors的个数
    for(i = n-1; i >= 0; --i) {
        int ifctmin = 0;      //ifctmin局部变量
        uint64_t min = fct[0]*factors[i];    //局部变量 最小为1*factors[i]
        for(ifct = 1; ifct < (int)nfct; ++ifct) {
            uint64_t k = fct[ifct]*factors[i];
            if(k < min) {
                min = k;
                ifctmin = ifct;  //ifctmin=ifct
            }
        }
        fct[ifctmin] *= factors[i];
    }
}
/** Returns the smallest integer bigger or equals to number that has a
    max prime factor smaller or equal than max_factor. If growth_limit
    is > 0, and no integer smaller than number+growth_limit is found,
    then the one with the smallest bigger prime factor between number
    and number+growth_limit is returned. */
uint64_t se_next_int_with_max_factor(uint64_t number, uint64_t max_factor, int64_t growth_limit)
{
    uint64_t result, candidate;
    uint64_t min = max_factor;  //max_factor=31

    for(result = candidate = number; 
        growth_limit <= 0 || candidate < number+growth_limit; ++candidate) {

        uint64_t factors[SE_MAX_N_FACTORS];
        int n_factors = se_factor(candidate, factors);

        if( factors[n_factors-1] <= max_factor ) {  //最大因子小于max_factor 返回最小的大于或者等于number的数，
//且这个数有一个小于或者等于max_factor的质因子 如果growth_limit>0 ，不会有比比number+growth_limit小的数
//被选取。然后，如果number太小 ok 如果 质因子>max_factor 


//返回挨着的 最小的整数 大于或者等于number 且这个整数有一个最大质因数 小于或者等于 max_factor
            result = candidate;
            break;
        }

        if( min > factors[n_factors-1] ) {
            min = factors[n_factors-1];
            result = candidate;
        }
    }

    return result;
}

