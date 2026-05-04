#include "../include/se_array_bit.h"


/* 2^ADDRESS_BITS_PER_UNIT == sizeof(bit_storage_unit) */
typedef int64_t bit_storage_unit;

#define ADDRESS_BITS_PER_UNIT 6
#define BITS_PER_UNIT (1<<ADDRESS_BITS_PER_UNIT)
#define BIT_INDEX_MASK (BITS_PER_UNIT - 1)

typedef struct bit_array_impl_s {
    /**
     * The bits in this bit_array.  The ith bit is stored in
     * bits[i/(8*sizeof(bit_storage_unit))] at bit position
     * i%(8*sizeof(bit_storage_unit)), where bit position 0 refers to the
     * least significant bit and sizeof(bit_storage_unit)-1 refers to
     * the most significant bit.
     */
    bit_storage_unit* bits;

    /** number of units allocated */
    size_t size;

    /** The number of units in the logical size of this bit_array.*/
    size_t units_in_use;

} bit_array_impl;

/**
 * Given a bit index return unit index containing it.
 */
#define unit_index(bitIndex) ((bitIndex) >> ADDRESS_BITS_PER_UNIT)
/*
static int unit_index(int bitIndex)
{
    int r1 = ((bitIndex) >> ADDRESS_BITS_PER_UNIT);
    int r2 = bitIndex / (8*sizeof(bit_storage_unit));
    if(r1 != r2) INFO(("unit_index(%d)=%d, %d\n", bitIndex, r1, r2));
    return r1;
}
*/


/**
 * Given a bit index, return a unit that masks that bit in its unit.
 */
#define bit(bitIndex) ((bit_storage_unit)(((bit_storage_unit)1) << ((bitIndex) & BIT_INDEX_MASK)))
/*
static bit_storage_unit bit(int bitIndex)
{
    bit_storage_unit r1 = ((((bit_storage_unit)1) << ((bitIndex) & BIT_INDEX_MASK)));
    bit_storage_unit r2 = 1<<(bitIndex%(8*sizeof(bit_storage_unit)));
    if(r1 != r2) INFO(("       bit(%d)=%lx, %lx\n", bitIndex, r1, r2));
    return r1;
}
*/

/**
 * Creates a bit set with initial size large enough to explicitly
 * represent bits with indices in the range 0 through nbits-1.
 * All bits are initially false.
 *
 * If nbits == 0, a default initial size will be used.
 */
bit_array ba_create(size_t nbits)
{
    bit_array_impl* ba = (bit_array_impl *)malloc(sizeof(*ba));
    memset(ba, 0, sizeof(*ba));

    if (nbits == 0) nbits = 4*BITS_PER_UNIT;

    ba->size = ( unit_index(nbits-1) + 1 );
    ba->bits = (bit_storage_unit *)malloc(  sizeof(bit_storage_unit) * ba->size );
    memset( ba->bits, 0, sizeof(bit_storage_unit) * ba->size );
    
    return ba;
}

/** Release used resources */
void ba_destroy(bit_array _ba)
{
    bit_array_impl* ba = (bit_array_impl*)_ba;
    free(ba->bits);
    free(ba);
}

byte* ba_serialize(bit_array _ba, int destroy, int64_t* nbytes)
{
    bit_array_impl* ba = (bit_array_impl*)_ba;
    int64_t size = sizeof(int64_t) + ba->units_in_use * sizeof(bit_storage_unit);
    int64_t inuse = ba->units_in_use;
    byte* buff;
    if(destroy) {
        buff = (byte *)realloc(ba->bits, (size_t)size);
        memmove(buff+sizeof(int64_t), buff, ba->units_in_use * sizeof(bit_storage_unit));
        memset(ba, 0, sizeof(bit_array_impl));/*safety precaution, in case someone will try to use it again*/
        free(ba);
    } else {
        buff = (byte *)malloc((size_t)size);
        memcpy(buff+sizeof(int64_t), ba->bits, ba->units_in_use * sizeof(bit_storage_unit));
    }

    memcpy(buff, &inuse, sizeof(int64_t));
    order_bytes_4(buff, (size_t)size, little_endian);

    *nbytes = size;
    return buff;
}

bit_array ba_deserialize(byte* buff, int allocmem)
{
    bit_array_impl* ba;
    int64_t inuse;
    
    ba = (bit_array_impl *)malloc(sizeof(*ba));

    memcpy(&inuse, buff, sizeof(int64_t));
    order_bytes_4((byte*)&inuse, sizeof(int64_t), little_endian);
    
    ba->units_in_use = ba->size = (size_t)inuse;
    if(allocmem) {
        ba->bits = (bit_storage_unit *)malloc(sizeof(bit_storage_unit) * ba->size);
        memcpy(ba->bits, buff+sizeof(int64_t), sizeof(bit_storage_unit) * ba->size);
    } else {
        memmove(buff, buff+sizeof(int64_t), ba->size * sizeof(bit_storage_unit));
        ba->bits = (bit_storage_unit*)buff;
    }
    order_bytes_4((byte*)ba->bits, ba->units_in_use * sizeof(bit_storage_unit), little_endian);
    return ba;
}

/**
 * Set the field units_in_use with the logical size in units of the
 * bit array.
 *
 * WARNING:This function assumes that the number of units actually in
 * use is less than or equal to the current value of units_in_use!
 */
static void ba_recalculate_units_in_use(bit_array_impl* ba)
{
    if(ba->units_in_use) {
        size_t i = ba->units_in_use;
        while(i > 0) {
            if(ba->bits[--i] != 0) {
                ba->units_in_use = i+1;
                return;
            }
        }
        ba->units_in_use = 0;
    }
}

/**
 * Ensures that the bit_array can hold enough units.
 */
static void ba_ensure_capacity(bit_array_impl* ba, int units_required)
{
    if((int)ba->size < units_required) {
        bit_storage_unit* newbits;
        size_t request = 2*ba->size;
        if((int)request < units_required) {
            request = units_required;
        }

        newbits = (bit_storage_unit *)realloc(ba->bits, sizeof(bit_storage_unit) * request);
        memset( newbits + ba->units_in_use, 0, 
                sizeof(bit_storage_unit) * (request - ba->units_in_use) );

        ba->bits = newbits;
        ba->size = request;
    }
}


/**
 * Sets the bit at the specified index to true.
 */
void ba_set(bit_array _ba, size_t bit_index)
{
    bit_array_impl* ba = (bit_array_impl*)_ba;
    size_t ui, ur;

    ui = unit_index(bit_index);
    ur = ui + 1;

    if(ba->units_in_use < ur) {
        ba_ensure_capacity(ba, ur);
        ba->units_in_use = ur;
    }
    ba->bits[ui] |= bit(bit_index);
}

/**
 * Sets the bit specified by the index to false.
 */
void ba_clear(bit_array _ba, size_t bit_index)
{
    bit_array_impl* ba = (bit_array_impl*)_ba;
    size_t ui;

    ui = unit_index(bit_index);
    if( ui >= ba->units_in_use )
        return;

    ba->bits[ui] &= ~bit(bit_index);

    if(ba->bits[ba->units_in_use-1] == 0) {
        ba_recalculate_units_in_use(ba);
    }
}

/**
 * Sets all the bits to false.
 */
void ba_clear_all(bit_array _ba)
{
    bit_array_impl* ba = (bit_array_impl*)_ba;
    while(ba->units_in_use > 0)
        ba->bits[--ba->units_in_use] = 0;
}

/**
 * Returns the value of the bit with the specified index. The value is
 * true if the bit with the index bit_index is currently set;
 * otherwise, the result is false.
 */
int ba_get(const bit_array _ba, size_t bit_index)
{
    const bit_array_impl* ba = (bit_array_impl*)_ba;
    size_t ui;

    ui = unit_index(bit_index);
    if( ui < ba->units_in_use )
        if( 0 != (ba->bits[ui] & bit(bit_index)) )
            return 1;

    return 0;
}

/**
 * Performs a logical OR of arguments. The result is stored in the
 * first parameter.
 */
void ba_or(bit_array _ba1, const bit_array _ba2)
{
    bit_array_impl* ba1 = (bit_array_impl*)_ba1;
    bit_array_impl* ba2 = (bit_array_impl*)_ba2;
    size_t units_in_common, i;

    if (ba1 == ba2)
        return;

    ba_ensure_capacity(ba1, ba2->units_in_use);

    /* Perform logical OR on bits in common */
    if(ba1->units_in_use < ba2->units_in_use) {
        units_in_common = ba1->units_in_use;
    } else {
        units_in_common = ba2->units_in_use;
    }

    for(i = 0; i < units_in_common; ++i)
        ba1->bits[i] |= ba2->bits[i];

    /* Copy any remaining bits */
    for(; i < ba2->units_in_use; ++i)
        ba1->bits[i] = ba2->bits[i];

    if (ba1->units_in_use < ba2->units_in_use)
        ba1->units_in_use = ba2->units_in_use;
}


static void ba_print_bits(FILE* f, bit_storage_unit b, int pidx)
{
    int i;
    int bits[8*sizeof(b)];

    for(i = 8*sizeof(b)-1; i >= 0; --i) {
        bit_storage_unit m = ((bit_storage_unit)1)<<i;
        if(0 != (m&b)) bits[i] = 1;
        else           bits[i] = 0;
        if(pidx) fprintf(f, "%2d ", i);
    }

    if(pidx) fprintf(f, "\n");

    for(i = 8*sizeof(b)-1; i >= 0; --i) {
        fprintf(f, "%2d ", bits[i]);
    }
    fprintf(f, "\n");
}

void ba_status(bit_array _ba)
{
    bit_array_impl* ba = (bit_array_impl*)_ba;
    int i;
    INFO(("There are %d / %d units in use:\n", ba->units_in_use, (int)ba->size));
    for(i=0; i<(int)ba->units_in_use; i++) {
        ba_print_bits(stderr, ba->bits[i], !i);
    }
}
