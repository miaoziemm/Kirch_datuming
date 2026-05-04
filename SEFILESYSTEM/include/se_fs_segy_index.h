#ifndef SE_FS_SEGY_INDEX_H
#define SE_FS_SEGY_INDEX_H

#include <SEBASIC/include/se_basic.h>
#include "se_fs_segy.h"
#include "se_fs_io.h"
#include "se_uri.h"

typedef void* segy_trace_index_t;

int segy_check_version_trace_index(const char* file, char** error);

segy_trace_index_t segy_open_trace_index(const char* file,
                                               const survey_description_t* sd,
                                               void* p_app_check,
                                               int (*trace_included)(void* p, double sx, double sy,
                                                                     double gx, double gy),
                                               void* p_app,
                                               int (*cache_index)(void* p,
                                                                  double x, double y),
                                               int sort_type);

void segy_close_trace_index(segy_trace_index_t _t);

int segy_trace_index_get_traces(segy_trace_index_t _t,
                                   double minof, double maxof,
                                   double mincx, double maxcx,
                                   double mincy, double maxcy,
                                   int* idx, int n);
void segy_index_sort_all_traces(void);
void segy_index_test_trace_index(void);


int segy_trace_index_get_traces_src( segy_trace_index_t _t,
                                        double minx, double maxx,
                                        double miny, double maxy,
                                        int** idx_ );

int segy_trace_index_get_traces_rec( segy_trace_index_t _t,
                                        double minx, double maxx,
                                        double miny, double maxy,
                                        int** idx_ );

int segy_trace_index_get_noffsets(segy_trace_index_t _t, double delta);

double segy_trace_index_get_delta_offset(segy_trace_index_t _t, int n);

int segy_trace_index_get_nazimuths(segy_trace_index_t _t, double delta);

double segy_trace_index_get_delta_azimuth(segy_trace_index_t _t, int n);

typedef struct segy_trace_index_statistics_s {
    double mino, maxo;
    double minox, maxox;
    double minoy, maxoy;
    double mina, maxa;
    double minsx, maxsx, mingx, maxgx, mincx, maxcx;
    double minsy, maxsy, mingy, maxgy, mincy, maxcy;
} segy_trace_index_statistics_t;

void segy_trace_index_getminmax(segy_trace_index_t _t,
                                   segy_trace_index_statistics_t* st);


typedef struct segy_trace_index_offset_s {
    double o;
    int64_t n;
} segy_trace_index_offset_t;

segy_trace_index_offset_t*
segy_trace_index_get_offset_distribution(segy_trace_index_t _t, int n);

typedef struct segy_trace_index_azimuth_s {
    double a;
    int64_t n;
} segy_trace_index_azimuth_t;

segy_trace_index_azimuth_t*
segy_trace_index_get_azimuth_distribution(segy_trace_index_t _t, int n);

typedef struct segy_trace_index_oa_s {
    double o;
    double a;
    int64_t n;
} segy_trace_index_oa_t;

segy_trace_index_oa_t*
segy_trace_index_get_oa_distribution(segy_trace_index_t _t, int no, int na);


typedef void* segy_trace_index_query_t;

segy_trace_index_query_t segy_trace_index_query_src(segy_trace_index_t _t,
                                                          int all,
                                                          double minx, double maxx,
                                                          double miny, double maxy);
segy_trace_index_query_t segy_trace_index_query_rec(segy_trace_index_t _t,
                                                          int all,
                                                          double minx, double maxx,
                                                          double miny, double maxy);
segy_trace_index_query_t segy_trace_index_query_cmp(segy_trace_index_t _t,
                                                          int all,
                                                          double minx, double maxx,
                                                          double miny, double maxy);

segy_trace_index_query_t
segy_trace_index_query_gather_src(segy_trace_index_t _t,
                                     double x, double y, double dx, double dy);

segy_trace_index_query_t
segy_trace_index_query_gather_rec(segy_trace_index_t _t,
                                     double x, double y, double dx, double dy);

segy_trace_index_query_t
segy_trace_index_query_gather_cmp(segy_trace_index_t _t,
                                     double x, double y, double dx, double dy);

void segy_trace_index_destroy_query(segy_trace_index_query_t _q);


typedef struct segy_trace_index_query_result_s {
    int64_t id;
    double x, y;
    double x1, y1;
} segy_trace_index_query_result_t;

/** Returns zero if there are no more results, in which case *res is left untouched. */
int segy_trace_index_query_next(segy_trace_index_query_t _q,
                                   segy_trace_index_query_result_t* res);


int segy_trace_index_get_close_src(segy_trace_index_t _t,
                                      double x, double y, double dx, double dy,
                                      segy_trace_index_query_result_t* res);

int segy_trace_index_get_close_rec(segy_trace_index_t _t,
                                      double x, double y, double dx, double dy,
                                      segy_trace_index_query_result_t* res);

int segy_trace_index_get_close_cmp(segy_trace_index_t _t,
                                      double x, double y, double dx, double dy,
                                      segy_trace_index_query_result_t* res);


/**
 * Returns the path to the trace index file; use this function
 * everywhere the index file is required, so the default value is
 * consistent everywhere.
 *
 * The returned value must be freed with free().
 */
char* segy_get_index_file(const char* segy_par_name);

/**
 * Returns the path to the trace coordinates file; use this function
 * everywhere the index file is required, so the default value is
 * consistent everywhere.
 *
 * The returned value must be freed with free().
 */
char* segy_get_coordinates_file(const char* segy_par_name);

void segy_create_trace_index(const char* file,
                                int is_coordinate_file,
                                const char* outfile,
                                const survey_description_t* sd,
                                const char* parset);

void segy_create_trace_index_v4(const char* file,
                                   int is_coordinate_file,
                                   const char* outfile,
                                   const survey_description_t* sd,
                                   const char* parset);

int segy_index_compute_offset_distribution(segy_trace_index_t _t,
                                             double min, double max, int nbins,
                                             double **values,
                                             int **counts,
                                             int* nvalues);


#endif