#include "../include/se_byte.h"

#if defined __GLIBC__ && ! defined(__INTEL_COMPILER) && ! defined(ppc64)
#include <byteswap.h>
#endif


void order_bytes_2(byte* buff, const size_t n, const byte_order bo_input)
{
    ASSERTM( n%2 == 0, "The size of the input array should be multiple of 2" );

#if defined(USE_BYTESWAP_REF_IMPL_FOR_UNALIGNED)
    if( ( ((size_t)buff) % 2 ) != 0 ) {
        reference_order_bytes_2( buff, n, bo_input );
        return;
    }
#endif

    if( NEEDS_SWAPPING(bo_input) ) {

        const size_t N = n>>1; /* the number of 16-bit integers = n/2 */
        uint16_t* buff16 = (uint16_t*)buff;
        size_t i;

#if defined(__INTEL_COMPILER) && defined(__i386)
        for(i = 0; i < N; i++) { /* TO_DO: find an optimized way */
            uint16_t x = buff16[i];
            buff16[i] = (uint16_t)bswap_const_16(x);
        }
/* GNU compiler have bswap_16 */
#elif defined __GLIBC__ && ! defined(__INTEL_COMPILER) && ! defined(ppc64)

        for(i = 0; i < N; i++) {
            buff16[i] = bswap_16(buff16[i]);
        }

#else /* Generic code for swapping */

        for(i = 0; i < N; i++) {
            uint16_t x = buff16[i];
            buff16[i] = (uint16_t)bswap_const_16(x);
        }
#endif

    } /* END if need to be swapped */
}



void order_bytes_4(byte* buff, const size_t n, const byte_order bo_input)
{
    ASSERTM( n%4 == 0, "The size of the input array should be multiple of 4" );

#if defined(USE_BYTESWAP_REF_IMPL_FOR_UNALIGNED)
    if( ( ((size_t)buff) % 4 ) != 0 ) {
        reference_order_bytes_4( buff, n, bo_input );
        return;
    }
#endif

    if( NEEDS_SWAPPING(bo_input) ) {

        const size_t N = n>>2; /* the number of 32-bit integers = n/4 */
        uint32_t* buff32 = (uint32_t*)buff;
        size_t i;

#if defined(__INTEL_COMPILER) && defined(__i386)
        for(i = 0; i < N; i++) { /* Intel x86 define intrinsic _bswap */
            buff32[i] = _bswap(buff32[i]);
        }

/* GNU compiler have bswap_32 */
#elif defined __GLIBC__ && ! defined(__INTEL_COMPILER) && ! defined(ppc64)

        for(i = 0; i < N; i++) {
            buff32[i] = bswap_32(buff32[i]);
        }

#else /* Generic code for swapping */

        for(i = 0; i < N; i++) {
            uint32_t x = buff32[i];
            buff32[i] = bswap_const_32(x);
        }
#endif

    } /* END if need to be swapped */
}


void order_bytes_8(byte* buff, const size_t n, const byte_order bo_input)
{
    ASSERTM( n%8 == 0, "The size of the input array should be multiple of 8" );

#if defined(USE_BYTESWAP_REF_IMPL_FOR_UNALIGNED)
    if( ( ((size_t)buff) % 8 ) != 0 ) {
        reference_order_bytes_8( buff, n, bo_input );
        return;
    }
#endif

    if( NEEDS_SWAPPING(bo_input) ) {

        const size_t N = n>>2; /* the number of 32-bit integers = n/4 */
        uint32_t* buff32 = (uint32_t*)buff;
        size_t i;

#if defined(__INTEL_COMPILER) && defined(__i386)
        for(i = 0; i < N; i += 2) { /* Intel x86 define intrinsic _bswap */
            uint32_t x0 = buff32[i+0];
            uint32_t x1 = buff32[i+1];
            buff32[i+0] = _bswap(x1);
            buff32[i+1] = _bswap(x0);
        }

/* GNU compiler have bswap_16 */
#elif defined __GLIBC__ && ! defined(__INTEL_COMPILER) && ! defined(ppc64)

        for(i = 0; i < N; i += 2) {
            uint32_t x0 = buff32[i+0];
            uint32_t x1 = buff32[i+1];
            buff32[i+0] = bswap_32(x1);
            buff32[i+1] = bswap_32(x0);
        }

#else /* Generic code for swapping */

        for(i = 0; i < N; i += 2) {
            uint32_t x0 = buff32[i+0];
            uint32_t x1 = buff32[i+1];
            buff32[i+0] = bswap_const_32(x1);
            buff32[i+1] = bswap_const_32(x0);
        }
#endif

    } /* END if need to be swapped */
}


void reference_order_bytes_2( byte* buff,
                                  const size_t n,
                                  const byte_order bo_input )
{
    ASSERTM( n%2 == 0, "The size of the input array should be multiple of 2" );

    if( NEEDS_SWAPPING(bo_input) ) {
        size_t i;
        for(i = 0; i < n; i += 2) {
            const size_t i0 = i + 0;
            const size_t i1 = i + 1;
            byte t;

            t = buff[i0]; buff[i0] = buff[i1]; buff[i1] = t;
        }
    }
}

void reference_order_bytes_4( byte* buff,
                                  const size_t n,
                                  const byte_order bo_input )
{
    ASSERTM( n%4 == 0, "The size of the input array should be multiple of 4" );

    if( NEEDS_SWAPPING(bo_input) ) {
        size_t i;
        for(i = 0; i < n; i += 4) {
            const size_t i0 = i + 0;
            const size_t i1 = i + 1;
            const size_t i2 = i + 2;
            const size_t i3 = i + 3;
            byte t;

            t = buff[i0]; buff[i0] = buff[i3]; buff[i3] = t;
            t = buff[i1]; buff[i1] = buff[i2]; buff[i2] = t;
        }
    }
}

void reference_order_bytes_8( byte* buff,
                                  const size_t n,
                                  const byte_order bo_input )
{
    ASSERTM( n%8 == 0, "The size of the input array should be multiple of 8" );

    if( NEEDS_SWAPPING(bo_input) ) {
        size_t i;
        for(i = 0; i < n; i += 8) {
            const size_t i0 = i + 0;
            const size_t i1 = i + 1;
            const size_t i2 = i + 2;
            const size_t i3 = i + 3;
            const size_t i4 = i + 4;
            const size_t i5 = i + 5;
            const size_t i6 = i + 6;
            const size_t i7 = i + 7;
            byte t;

            t = buff[i0]; buff[i0] = buff[i7]; buff[i7] = t;
            t = buff[i1]; buff[i1] = buff[i6]; buff[i6] = t;
            t = buff[i2]; buff[i2] = buff[i5]; buff[i5] = t;
            t = buff[i3]; buff[i3] = buff[i4]; buff[i4] = t;
        }
    }
}


void ibm_to_ieee_sep(const unsigned int* from, 
                         unsigned int* to, int n, int endian)
{
    int fconv, fmant, i, t, watchdog;

    for (i=0;i<n;++i) {

        fconv = from[i];

        /* if little endian, i.e. endian=0 do this */
        if (endian == 0)
            fconv = (fconv << 24) |
                ((fconv>>24)&0xff) |
                ((fconv&0xff00)<<8) |
                ((fconv&0xff0000)>>8);

        if (fconv) {
            fmant = 0x00ffffff & fconv;
            t = (int) ((0x7f000000 & fconv) >> 22) - 130;
            watchdog = 70;
            while (--watchdog && !(fmant & 0x00800000)) {
                --t;
                fmant <<= 1;
            }
            if( watchdog > 0) {
                if (t > 254) fconv = (0x80000000 & fconv) | 0x7f7fffff;
                else if (t <= 0) fconv = 0;
                else fconv = (0x80000000 & fconv) |
                         (t << 23)|(0x007fffff & fmant);
            } else {
                fconv = 0x7fc00000; /* NaN */
            }
        }
        to[i] = fconv;
    }
    return;
}

void ieee_to_ibm_sep(const unsigned int* from,
                         unsigned int* to, int n, int endian)
{
    int fconv, fmant, i, t, watchdog;

    for(i = 0; i < n; i++) {
	fconv = from[i];
	if( fconv ) {
	    fmant = (0x007fffff & fconv) | 0x00800000;
	    t = (int) ((0x7f800000 & fconv) >> 23) - 126;
            watchdog = 70;
	    while(--watchdog && (t & 0x3) ) {
                ++t;
                fmant >>= 1;
            }

            if( watchdog > 0) {
                fconv = (0x80000000 & fconv) | (((t>>2) + 64) << 24) | fmant;
            } else {
                fconv = 0x7fc00000; /* NaN (not sure for IBM) */
            }
	}

	if( endian == 0 ) {
            fconv = (fconv<<24) |
                ((fconv>>24)&0xff) |
                ((fconv&0xff00)<<8) |
                ((fconv&0xff0000)>>8);
        }

	to[i] = fconv;
    }
}




/* new, untested code  */


/**
   description :
      routine to convert ibm-floating point to
      machine-internal representation

   implementation :

     ibm-format

     seeeeeeemmmmmmmmmmmmmmmmmmmmmmmm

     s -sign
     e -exponent
     m -mantissa

     value of number is mantissa *16**(exponent-63)
     radix is to left of most significant mantissa digit

   input arguments :
     iar  -  input array
     n    -  number of elements

   output arguments :
     kar  -  output array

   local variables :
   problems :
     iar and far may have same start-address
*/
void ibm_to_ieee(const unsigned int *iar, unsigned int *kar, 
                     int n, int endian)
{
    int in;
    unsigned ic,ipos,iexp,man,ish;
    unsigned mask1,mask2,mask3,mask4,mask5;

    mask1 = 0x80000000;
    mask2 = 0x7fffffff;
    mask3 = 0x00ffffff;
    mask4 = 0x00800000;
    mask5 = 0x007fffff;

    for(in = 0; in < n; in++) {
        /* copy input into ic */
        ic = iar[in];
        if(endian == 0)
            ic = (ic << 24) | ((ic>>24)&0xff) | 
                ((ic&0xff00)<<8) | ((ic&0xff0000)>>8);

        /* put sign bit into output */
        kar[in] = ic & mask1;
        /* mask off sign bit put rest into ipos */
        ipos = ic & mask2;
        /* get exponent into iexp */
        iexp = ipos >> 24;
        /* mask off mantissa and put into man */
        man = ipos & mask3;
        /* set a counter to 1 */
        ish = 1;
        /* start shifting the mantissa */
        if((man & mask4) == 0) { man = man << 1; ish++; }
        if((man & mask4) == 0) { man = man << 1; ish++; }
        if((man & mask4) == 0) { man = man << 1; ish++; }
        /* add in excess 127 and shift exponent to proper place */
        iexp = ((((iexp - 64) << 2) - ish + 127) << 23);
        /* or the pieces together */
        kar[in] = ((man & mask5) | (kar[in] | iexp));
        if(ic == 0) kar[in] = 0;
    }
}


/**
   description :
     routine to convert ieee-floating point (sun) to
     ibm floating point

   implementation :

     ibm-format

     seeeeeeemmmmmmmmmmmmmmmmmmmmmmmm

     s -sign
     e -exponent (7 bit)
     m -mantissa (24 bit)

     value of number is mantissa *16**(exponent-64)
     radix is to left of most significant mantissa digit

     ieee-format

     seeeeeeeemmmmmmmmmmmmmmmmmmmmmmm

     s -sign
     e -exponent (8 bit)
     m -mantissa (23 bit)

     value of number is mantissa *2**(exponent-127)
     the 23 bit mantissa has an implicit 1 bit in the most
     significant position to make it 24 bits. radix after
     implicit bit

   referenced subroutines :

   input arguments :
     iar  -  input array
     n    -  number of elements

  output arguments :
     kar  -  output array

   local variables :

   problems :
     iar and kar may have same start-address
*/
void ieee_to_ibm(const unsigned *iar, unsigned int *kar, int n, int endian)
{
    int in;
    unsigned ic,oc,isign,icc,iover,xp,l;
    unsigned z8,z008;
    unsigned mask,mask2,mask3,mask4,mask5;

    z8 =    0x80000000;
    z008 =  0x00800000;

    mask =  0x00000003;
    mask2 = 0x00ffffff;
    mask3 = 0xff000000;
    mask4 = 0x000000ff;
    mask5 = 0x7fffffff;

    for( in = 0; in < n; in++) {
        /* copy input into ic */
        ic = iar[in];
        /* sign bit of numbers into isign */
        isign = ic & z8;
        /* shift right and form exponent */
        icc = (((ic & mask5) >> 23) & mask4)+133;
        /* add in implied bit */
        l = (mask2 & (ic | z008));
        /* now adjust exponent to base 16 and also the mantissa. */
        /* this is done by shifting mantissa mantissa exp-exp/4*4 times. */
        iover = (icc & mask);
        if((3-(int)iover) > 0) l = l >> (3-(int)iover);
        if((3-(int)iover) < 0) l = l << (3-(int)iover);
        xp = ((icc << 22) & mask3);
        /* or the pieces together */
        oc = (l | (xp | isign));
        if(ic == 0) oc = 0;
	if( endian == 0 ) 
            oc = (oc<<24) | ((oc>>24)&0xff) | 
                ((oc&0xff00)<<8) | ((oc&0xff0000)>>8);
        kar[in] = oc;
    }
}
