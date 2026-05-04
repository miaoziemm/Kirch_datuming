#ifndef SE_FS_SEGY_H
#define SE_FS_SEGY_H
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include "se_fs_segy_header.h"
#include "se_fs_segy_bin_types.h"
#include "se_fs_segy_index.h"
#include "se_fs_segy_trchead_eval.h"
#include "se_fs_io.h"
#include "se_fs_sep.h"

typedef void* segy_read_t;

typedef struct segy_reader_params_s {
    int    use_scalco;
    double scalco_factor;

    int ns;
    int forced_format;
    int forced_ns;
    int forced_idt;

    char* sx_str;
    char* sy_str;
    char* gx_str;
    char* gy_str;

    trchead_eval_handle sx;
    trchead_eval_handle sy;
    trchead_eval_handle gx;
    trchead_eval_handle gy;

    char* weights;
    sep_t* sepweights;

    int use_cache;
    size_t hashsize;
    int64_t max_cache_size;

    int keep_open;

    int use_mmap;
    char* mmap_flags;
} segy_reader_params_t;


inline void segy_init_reader_par(segy_reader_params_t* p)
{
    memset(p, 0, sizeof(*p));
    p->use_scalco = 1;
    p->scalco_factor = 1.0;
    p->use_cache = 0;
    p->hashsize = 4096;
    p->max_cache_size = 20*1024*1024;

    p->keep_open = 0;
    p->use_mmap = 0;
}


segy_read_t segy_create_reader( const char* file_name,
                                      segy_reader_params_t* p);
segy_read_t segy_get_segy_info( const char* file_name,
                                      int* ns, int* dt,
                                      int* ntraces,
                                      const char* parset);

void segy_init_reader_par_from_parset(segy_reader_params_t* p, const char* parset);
void segy_init_reader_par_expression_strings(segy_reader_params_t* p,
                                                const char* sx,
                                                const char* sy,
                                                const char* gx,
                                                const char* gy);
void segy_init_reader_par_expressions(segy_reader_params_t* p);
void segy_free_reader_par(segy_reader_params_t* p);


int segy_get_segy_bps(const segy_read_t segy );
int segy_get_trace_length(const segy_read_t segy );
const byte* segy_get_reel_header(const segy_read_t segy);
char* segy_get_ebcdic_text(const segy_read_t segy);
void segy_get_bhead(const segy_read_t segy, binhead* bhead);

int segy_get_ns(const segy_read_t segy);
int segy_get_ntraces(const segy_read_t segy);
int segy_get_dt(const segy_read_t segy);
int segy_get_format(const segy_read_t segy);
int segy_get_nfiles(const segy_read_t segy);

void segy_destroy_segy_reader(segy_read_t segy);


void segy_compute_trace_coordinates( const segy_read_t segy,
                                        const trchead* head, 
                                        float* lsx, float* lsy,
                                        float* lgx, float* lgy );


void segy_coords_to_integer(const segy_read_t segy,
                               const trchead* head,
                               double sx, double sy, double gx, double gy,
                               int* isx, int* isy, int* igx, int* igy);

void segy_ieee2segy(const segy_read_t segy,
                       trchead* head,
                       byte* bdata, const float* data,
                       int nsamples);


int segy_get_segy_traces( segy_read_t segy,
                             int first_trace, int ntraces, const int* index_traces, 
                             float* data, float* coords,  int coords_only);
int segy_get_raw_segy_traces( segy_read_t segy,
                                 int first_trace, int ntraces,
                                 const int* index_traces,
                                 byte* data);

int segy_get_segy_traces_and_headers( segy_read_t segy,
                                         int first_trace, 
                                         int ntraces, const int* index_traces,
                                         float* data, float* coords, 
                                         trchead* hdrs );



void segy_binheader_order_bytes(binhead* HEADER);
void segy_trcheader_order_bytes(trchead* HEADER);



typedef void (*segy_fix_data_fnc)(const trchead* head, 
                                  const byte* bdata,
                                  float* data, int nsamples);
typedef void (*segy_ieee2segy_fnc)(trchead* head, 
                                   byte* bdata, const float* data,
                                   int nsamples);


void segy_fix_data_ibm_float( const trchead* head, 
                              const byte* bdata, float* data, int nsamples );
void segy_fix_data_ieee_float( const trchead* head, 
                               const byte* bdata, float* data, int nsamples );
void segy_fix_data_fixed_point_4B( const trchead* head, 
                                   const byte* bdata, float* data, 
                                   int nsamples );
void segy_fix_data_fixed_point_2B( const trchead* head, 
                                   const byte* bdata, float* data, 
                                   int nsamples );
void segy_fix_data_fixed_point_1B( const trchead* head, 
                                   const byte* bdata, float* data, 
                                   int nsamples );
void segy_fix_data_fixed_point_gain( const trchead* head, 
                                     const byte* bdata, float* data, 
                                     int nsamples );

void segy_ieee2segy_ibm_float( trchead* head, 
                               byte* bdata, const float* data, int nsamples );
void segy_ieee2segy_ieee_float( trchead* head, 
                                byte* bdata, const float* data, int nsamples );
void segy_ieee2segy_fixed_point_4B( trchead* head, 
                                    byte* bdata, const float* data, 
                                    int nsamples );
void segy_ieee2segy_fixed_point_2B( trchead* head, 
                                    byte* bdata, const float* data, 
                                    int nsamples );
void segy_ieee2segy_fixed_point_1B( trchead* head, 
                                    byte* bdata, const float* data, 
                                    int nsamples );
void segy_ieee2segy_fixed_point_gain( trchead* head, 
                                      byte* bdata, const float* data, 
                                      int nsamples );



inline double segy_compute_azimuth(double sx, double sy, double gx, double gy)
{
    /*atan2(sx-gx, sy-gy)/M_PI*180.0 + 180.0; ?? */
    double a = atan2(gx-sx,gy-sy);
    if(a < 0) a += 2*M_PI;
    return a*(180.0/M_PI);
}

inline double segy_compute_azimuth_coords(float* coords)
{
    return segy_compute_azimuth(coords[0], coords[1], coords[2], coords[3]);
}

int segy_detect_null_trace(const float* tr, int64_t n);

/**
 * Unified way of dealing with first_trace,last_trace,ntraces
 * parameters.  Returns NULL if successfull or an error message that
 * you need to free in case of error. The ouptut parameters
 * first_trace and last_trace are updated only in case of success.
 *
 * The values for first_trace and last_trace will be zero-based, even
 * if the user sees them as starting from one.
 *
 * To loop over the traces one can then use something like:
 * for(it = first_trace; it <= last_trace; ++it)
 */
char* segy_get_first_last_trace_parameters(int ntraces,
                                              int* first_trace, int* last_trace);


void segy_bin_check_and_convert_friendly_pars(void);



#endif