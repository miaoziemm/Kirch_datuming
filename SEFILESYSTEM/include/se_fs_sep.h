#ifndef SE_FS_SEP_H
#define SE_FS_SEP_H
#include <SEBASIC/include/se_basic.h>
#include "se_fs_io.h"
#include "se_uri.h"
#include "se_par_sep.h"

#define SEP_DEFAULT_LABEL_X "crosslines"
#define SEP_DEFAULT_LABEL_Y "inlines"

typedef struct sep_cube_s {
    int  n1, n2, n3, n4, n5, n6, n7, n8, n9;
    real d1, d2, d3, d4, d5, d6, d7, d8, d9;
    real o1, o2, o3, o4, o5, o6, o7, o8, o9;
} sep_cube_t;

#define SEP_NDIM_MAX             9    /* Max number of dimension in SEP */  //SEP中的最大维数
#define SEP_EOL                 12    /* ASCII EOL (Ctrl+L) */   //
#define SEP_EOT                  4    /* ASCII EOT (End of Transmission) */
#define SEP_MAX_NONTEXT_CHARS   12    /* non-text sep gheader char limit */

#define SEP_READ   0x01
#define SEP_WRITE  0x02

/* SEP global headers */
typedef struct sep_headers_s {
    char*    filename;          /* .H filename */  //文件名字
    int      fd;                /* .H file descriptor */  //文件描述符号 主题词
    int      ndim;              /* number of dimensions: 1 to SEP_NDIM_MAX */    //维数
    int      n[SEP_NDIM_MAX];/* size on each axis */   //坐标大小  //SEP_NDIM_MAX  最大维数
    double   d[SEP_NDIM_MAX];/* distance on each axis */
    double   o[SEP_NDIM_MAX];/* origin of each axis */
    char*    in;                /* .H@ filename */
    char*    expanded_in;       /* .H@ filename, expanded relative to the header */
    int      esize;             /* value size in bytes */
    int      le;                /* little endian flag */  //是否是little endian
    se_hash ht_int;            /* int global headers (name/int se_hash) */
    se_hash ht_float;          /* float global headers (name/float se_hash) */
    se_hash ht_str;            /* string global headers  (name/char* se_hash) */
    off_t    raw_size;          /* size (in bytes) of the headers */  //字节 头的字节
    char*    history;           /* old .H headers (in full) */
    int      has_inline_data;   /* the .H file (or stream) has inline data */
    int      is_modified;       /* is altered (since last write) flag */
    int      is_commited;       /* is written on file flag */
} sep_headers_t;

/* SEP data */
typedef struct sep_data_s {
    char*    filename;          /* .H@ filename */        
    se_fsio*  io;               /* i/o handler */
    int      ntraces;           /* trace count */
    int      trace_size;        /* number of values per trace */
    int      crt_trace;         /* current trace index */
    int      trace_index;       /* trace index */
    int      can_seek;          /* stream can seek flag */
    int      is_commited;       /* at least one trace was written flag */
    off_t    data_offset;       /* data offset, usualy zero */
    off_t    header_size;       /* trace header size, usualy zero */
} sep_data_t;


/* A SEP datasource */
typedef struct sep_s {
    sep_headers_t* headers;     /* SEP global headers (.H) */
    sep_data_t*    data;        /* SEP data (.H@) */
    int               mode;        /* access mode (eg. SEP_WRITE) */
    int               data_only;   /* no global headers (data only) */
    int               headers_only;
} sep_t;


sep_t* sep_open( const char* path_uri, int mode, int headers_only );

/** At the moment this is not 100% consistent with sep_open */
int sep_is_good_sep_header(const char* file);

void sep_close( sep_t* sep );
void sep_write_headers( sep_t* sep );
void sep_clear_history( sep_t* sep );
float* sep_read_fromsep(sep_t* sep, int64_t size, int64_t offset, float** buff, size_t* crt_size);
float* sep_read_fromfile(const char* file, int64_t size, int64_t offset);
float* sep_read_whole_file(const char* file);
float* sep_read_whole_file_check_values(const char* file, 
                                           float min, float max, float replace,
                                           int exit_on_error);

void sep_copy_headers(sep_t* dst, const sep_t* src);

void sep_fill_headers_from_par(sep_t* sep, const char* pargroup);

rsf3d_t* init_rsf3d_from_sep(const char* file);
rsf2d_t* init_rsf2d_from_sep(const char* file);

// grid for rsf need not match sep axes; returns 1 if sep grid does not cover rsf
int interp_rsf3d_from_sep(sep_t *sep, rsf3d_t *rsf);

// these two functions read the 3d cube starting at the current offset
float* sep_read_part_sep(const sep_t* sep,
                            const int m1, const int m2, const int m3,
                            const int k1, const int k2, const int k3);
rsf3d_t* init_rsf3d_part_from_sep(const sep_t* sep,
                                          const int m1, const int m2, const int m3,
                                          const int k1, const int k2, const int k3);

// these two functions read the 3d cube from specified offset off
float* sep_read_part_sep_offset(const sep_t* sep, off_t off,
                                   const int m1, const int m2, const int m3,
                                   const int k1, const int k2, const int k3);
rsf3d_t* init_rsf3d_part_from_sep_offset(const sep_t* sep, off_t off,
                                                 const int m1, const int m2, const int m3,
                                                 const int k1, const int k2, const int k3);

int sep_have_hdr_int(const sep_t* sep, const char* name);
int sep_get_hdr_int(const sep_t* sep, const char* name, int def);
int sep_have_hdr_float(const sep_t* sep, const char* name);
double sep_get_hdr_float(const sep_t* sep, const char* name, double def);
int sep_have_hdr(const sep_t* sep, const char* name);
char* sep_get_hdr(const sep_t* sep, const char* name, const char* def);

void sep_set_header(sep_t* sep, const char* name, const char* value);
void sep_set_header_int(sep_t* sep, const char* name, int value);
void sep_set_header_float(sep_t* sep, const char* name, double value);

void sep_set_header_survey_description(sep_t* sep,
                                          const survey_description_t* sd,
                                          const char* prefix,
                                          int use_only_set_fields);
void sep_get_header_survey_description(const sep_t* sep,
                                          survey_description_t* sd,
                                          const char* prefix,
                                          const survey_description_t* project_sd,
                                          int adjust_values);

void sep_set_axis(sep_t* sep, int aindex, int n, double o, double d, const char* label);
void sep_insert_axis(sep_t* sep, int aindex, int n, double o, double d, const char* label);

inline int64_t sep_get_total_size(const sep_t* sep)
{
    int64_t i, n = 1;
    for(i = 0; i < sep->headers->ndim; ++i) {
        n*= sep->headers->n[i];
     
    }
    return n;
}

inline int64_t sep_get_min_ndim(const sep_t* sep)
{
    int64_t i;
    for(i = sep->headers->ndim-1; i >= 0; --i) {
        if ( (sep->headers->n[i] > 1) || (!FISZERO(sep->headers->o[i])) ) {
            return i+1;
        }
    }
    return 1;
}

inline int sep_same_dimensions(const sep_t* sep1,
                                      const sep_t* sep2)
{
    int64_t i;
    if(sep1->headers->ndim > sep2->headers->ndim) {
        const sep_t* tmp = sep1;
        sep1 = sep2;
        sep2 = tmp;
    }

    for(i = 0; i < sep1->headers->ndim; ++i) {
        if(sep1->headers->n[i] != sep2->headers->n[i])
            return 0;
        if(!FEQUAL(sep1->headers->o[i], sep2->headers->o[i]))
            return 0;
        if(!FEQUAL(sep1->headers->d[i], sep2->headers->d[i]))
            return 0;
    }

    for(i = sep1->headers->ndim; i < sep2->headers->ndim; ++i) {
        if(sep2->headers->n[i] != 1)
            return 0;
    }

    return 1;
}

inline grid3d_t* init_grid3d_part_from_sepheaders(const sep_t* sep,
		const int m1, const int m2, const int m3,
		const int k1, const int k2, const int k3)
{
    return create_grid3d
    		(m2, sep->headers->d[1], sep->headers->o[1] + k2 * sep->headers->d[1],
    		 m3, sep->headers->d[2], sep->headers->o[2] + k3 * sep->headers->d[2],
    		 m1, sep->headers->d[0], sep->headers->o[0] + k1 * sep->headers->d[0]);
}

inline grid3d_t* init_grid3d_from_sepheaders(const sep_t* sep)
{
    return create_grid3d
        ( sep->headers->n[1], sep->headers->d[1], sep->headers->o[1],
          sep->headers->n[2], sep->headers->d[2], sep->headers->o[2],
          sep->headers->n[0], sep->headers->d[0], sep->headers->o[0] );
}

inline grid2d_t* init_grid2d_from_sepheaders(const sep_t* sep)
{
    return create_grid2d
        ( sep->headers->n[0], sep->headers->d[0], sep->headers->o[0],
          sep->headers->n[1], sep->headers->d[1], sep->headers->o[1] );
}

inline void init_grid_from_sepheaders(const sep_t* sep, int ndim, grid_t *grid)
{
    int i;

    ASSERT(ndim<=MAX_GRID_DIMMENSIONS);
    ASSERT(ndim<=SEP_NDIM_MAX);

    for(i=0;i<ndim;i++) {
        grid->a[i].n=sep->headers->n[i];
        grid->a[i].d=sep->headers->d[i];
        grid->a[i].o=sep->headers->o[i];
    }
    grid->ndim=ndim;
}

inline grid_t* create_grid_from_sepheaders(const sep_t* sep, int ndim)
{
    grid_t* grid = (grid_t*)malloc(sizeof(*grid));
    init_grid_from_sepheaders(sep, ndim, grid);
    return grid;
}

int sep_is_x_axis(const char* label);

int sep_is_y_axis(const char* label);

void sep_adjust_survey_description(const sep_t* sep,
                                       survey_description_t* sd);

inline void sep_get_indices(const sep_t* sep,  /* the SEP object */
                                    int64_t idx,          /* the linear index */
                                    int first_axis,       /* index of the first axis to consider */
                                    int last_axis,        /* index of the last axis to consider */
                                    int64_t* indices)     /* where to store the result: room for at least last_axis-first_axis+1 */
{
    const int naxis = last_axis - first_axis + 1;
    int iaxis;

    ASSERT(idx >= 0);
    ASSERT(first_axis > 0 && first_axis < SEP_NDIM_MAX);
    ASSERT(last_axis > 0 && last_axis < SEP_NDIM_MAX);
    ASSERT(first_axis <= last_axis);

    for(iaxis = 0; iaxis < naxis - 1; ++iaxis) {
        const int64_t n = sep->headers->n[first_axis + iaxis];
        indices[iaxis] = idx % n;
        idx /= n;
    }
    indices[naxis-1] = idx;
}



#endif