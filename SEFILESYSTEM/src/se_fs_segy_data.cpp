#include "../include/se_fs_segy.h"



void segy_binheader_order_bytes(binhead* HEADER)
{
    byte* h = (byte*)HEADER;

    order_bytes_4(h, 12, big_endian);
    order_bytes_2(h+12, 46, big_endian);
}


void segy_trcheader_order_bytes(trchead* HEADER)
{
    byte* h = (byte*)HEADER;
    size_t off, size;

    off = 0;     size = 28;
    order_bytes_4(h + off, size, big_endian);
    off += size; size = 8;
    order_bytes_2(h + off, size, big_endian);
    off += size; size = 32;
    order_bytes_4(h + off, size, big_endian);
    off += size; size = 4;
    order_bytes_2(h + off, size, big_endian);
    off += size; size = 16;
    order_bytes_4(h + off, size, big_endian);
    off += size; size = 92;
    order_bytes_2(h + off, size, big_endian);
    off += size; size = 28;
    order_bytes_4(h + off, size, big_endian);
}

void segy_fix_data_ibm_float( const trchead* head, 
                              const byte* bdata, float* data, int nsamples )
{
    (void)head;
    if(DEFAULT_BYTE_ORDER == little_endian) {
        ibm_to_ieee_sep((const unsigned int*)bdata,
                            (unsigned int*)data, nsamples, 0);
    } else {
        ibm_to_ieee_sep((const unsigned int*)bdata,
                            (unsigned int*)data, nsamples, 1);
    }
}

void segy_fix_data_ieee_float( const trchead* head, 
                               const byte* bdata, float* data, int nsamples )
{
    (void)head;
    memcpy(data, bdata, nsamples<<2);
    order_bytes_4((byte*)data, nsamples<<2, big_endian);
}


void segy_fix_data_fixed_point_4B( const trchead* head, 
                                   const byte* bdata, float* data, 
                                   int nsamples )
{
    int32_t* input = (int32_t*)data;
    double factor;

    if( head->trwf > 0 )
        factor = (double)1.0 / (double)( 0x0001 << (head->trwf) );
    else
        factor = (double)( 0x0001 << ( - head->trwf) );

    memcpy(data, bdata, nsamples<<2);
    order_bytes_4((byte*)data, nsamples<<2, big_endian);

    while(nsamples--) {
        int32_t i = *(input++);
        *(data++) = (float)(i*factor);
    }
}

void segy_fix_data_fixed_point_2B( const trchead* head, 
                                   const byte* bdata, float* data, 
                                   int nsamples )
{
    byte* bcopy = (byte*)(data + nsamples/2);
    double factor;

    if( head->trwf > 0 )
        factor = (double)1.0 / (double)( 0x0001 << (head->trwf) );
    else
        factor = (double)( 0x0001 << ( - head->trwf) );

    memcpy(bcopy, bdata, nsamples<<1);
    order_bytes_2(bcopy, nsamples<<1, big_endian);

#if defined(NEEDS_MEMALIGN)
    if( ( ((size_t)bcopy) %2 ) != 0 ) {
        union {
            byte b[2];
            int16_t i;
        }i;
        while(nsamples--) {
            i.b[0] = *(bcopy++);
            i.b[1] = *(bcopy++);
            *(data++) = (float)( i.i * factor);
        }
    } else 
#endif
    {
        int16_t* input = ((int16_t*)bcopy);
        while(nsamples--) {
            int16_t i = *(input++);
            *(data++) = (float)(i*factor);
        }
    }
}

void segy_fix_data_fixed_point_1B( const trchead* head, 
                                   const byte* bdata, float* data, 
                                   int nsamples )
{
    int8_t* input = (int8_t*)bdata;
    double factor;

    if( head->trwf > 0 )
        factor = (double)1.0 / (double)( 0x0001 << (head->trwf) );
    else
        factor = (double)( 0x0001 << ( - head->trwf) );

    while(nsamples--) {
        int8_t i = *(input++);
        *(data++) = (float)(i*factor);
    }
}

void segy_fix_data_fixed_point_gain( const trchead* head, 
                                     const byte* bdata, float* data, 
                                     int nsamples )
{
    (void)head;
    (void)bdata;
    (void)data;
    (void)nsamples;
    LOG_CONSOLE(SE_ERROR, "Data in 4B fixed point gain not supported yet");

}


void segy_ieee2segy_ibm_float( trchead* head, 
                               byte* bdata, const float* data, int nsamples )
{
    (void)head;

    if(DEFAULT_BYTE_ORDER == little_endian) {
        ieee_to_ibm((const unsigned int*)data,
                        (unsigned int*)bdata, nsamples, 0);
    } else {
        ieee_to_ibm((const unsigned int*)data,
                        (unsigned int*)bdata, nsamples, 1);
    }
}

void segy_ieee2segy_ieee_float( trchead* head, 
                                byte* bdata, const float* data, int nsamples )
{
    (void)head;
    memcpy(bdata, data, nsamples<<2);
    order_bytes_4(bdata, nsamples<<2, big_endian);
}

void segy_ieee2segy_fixed_point_4B( trchead* head, 
                                    byte* bdata, const float* data, 
                                    int nsamples )
{
    (void)head;
    (void)bdata;
    (void)data;
    (void)nsamples;
    LOG_CONSOLE(SE_ERROR, "Transforming data to 4B fixed point not supported yet");
}

void segy_ieee2segy_fixed_point_2B( trchead* head, 
                                    byte* bdata, const float* data, 
                                    int nsamples )
{
    (void)head;
    (void)bdata;
    (void)data;
    (void)nsamples;
    LOG_CONSOLE(SE_ERROR, "Transforming data to 2B fixed point not supported yet");
}

void segy_ieee2segy_fixed_point_1B( trchead* head, 
                                    byte* bdata, const float* data, 
                                    int nsamples )
{
    (void)head;
    (void)bdata;
    (void)data;
    (void)nsamples;
    LOG_CONSOLE(SE_ERROR, "Transforming data to 1B fixed point not supported yet");
}

void segy_ieee2segy_fixed_point_gain( trchead* head, 
                                      byte* bdata, const float* data, 
                                      int nsamples )
{
    (void)head;
    (void)bdata;
    (void)data;
    (void)nsamples;
    LOG_CONSOLE(SE_ERROR, "Transforming data to 4B fixed point gain not supported yet");
}
