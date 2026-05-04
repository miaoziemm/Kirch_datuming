#include "../include/se_fs_segy_index.h"
#include "../include/se_fs_sep.h"
#include "../include/se_jcs.h"

#include <stdlib.h>
#include <string.h>


#define SEGYTI_SCHEMA_VERSION 3
/*#define SEGYTI_SEPARATE_OFFSET_INDEX */
/* #define SEGYTI_USE_AZIMUTH_OFFSET_INDEX */

typedef struct segy_index_offset_bin_s {
    float off;
    int cnt;
}segy_index_offset_bin_t;

/* this structure is used both when creating a index file and when
   using it; it was just convenient to reuse some of the utility
   functions.  The meaning of few fields change, but this is handled
   in the corresponding routines.*/
typedef struct segyti_trace_index_s {
    char* file;
    char* parset;
    char* dbfile;
    int is_coordinate_file;
    survey_description_t sd;

    segy_read_t segy;
    sep_t* sep;

    int ntraces;

    int sort_type;
    int query_idx_minof;
    int query_idx_maxof;
    int query_idx_mincx;
    int query_idx_maxcx;
    int query_idx_mincy;
    int query_idx_maxcy;

    int nbuf;
    float* coords;

    sqlite3 *db;
    sqlite3_stmt* stmt_update_stat;
    sqlite3_stmt* stmt_update_minmax;
    sqlite3_stmt* stmt_update_xy;
    sqlite3_stmt* stmt_get_traces;

    int first_idx_in_db;
    segy_trace_index_statistics_t mm;

    timer* timer_reading;
    timer* timer_updating;
    timer* timer_optimization;

    void* p_app_check;
    int (*trace_included)(void* p, double sx, double sy, double gx, double gy);

    void* p_app;
    int (*cache_index)(void* p, double x, double y);
    int has_crt_data, get_traces_binded;
    double last_minof, last_maxof;
    double last_mincx, last_maxcx;
    double last_mincy, last_maxcy;


    /* Offset distribution */
    size_t offsets_bs;
    size_t noffsets;
    segy_index_offset_bin_t* offsets;
}segyti_trace_index_t;


#define DO_SQLITE(work) do {                                            \
        int rc = work;                                                  \
        if( rc != SQLITE_OK ) {                                         \
            char* msg = strdup(sqlite3_errmsg(t->db));              \
            sqlite3_close(t->db);                                       \
            ERROR(("Problems working with index file %s - %s returned: %s\n", \
                   t->dbfile, #work, msg));                             \
            free(msg);                                              \
        }                                                               \
    } while(0)

#define DO_INSERT(work) do {                                            \
        int rc = work;                                                  \
        if(rc != SQLITE_DONE) {                                         \
            char* msg = strdup(sqlite3_errmsg(t->db));              \
            sqlite3_close(t->db);                                       \
            ERROR(("Problems working with index file %s - %s returned: %s\n", \
                   t->dbfile, #work, msg));                             \
            free(msg);                                              \
        }                                                               \
    } while(0)

static void segyti_dosql(segyti_trace_index_t* t, const char* stmt, ...)
{
    va_list args;
    char* zSQL;

    va_start( args, stmt );
    zSQL = sqlite3_vmprintf(stmt, args);
    if( sqlite3_exec(t->db, zSQL, NULL, NULL, NULL) != SQLITE_OK ) {
        char* msg = strdup(sqlite3_errmsg(t->db));
        sqlite3_close(t->db);
        ERROR(("Problems working with index file %s - %s returned: %s",
               t->dbfile, zSQL, msg));
        free(msg);
    }
    sqlite3_free(zSQL);
}

static void segyti_update_stat(segyti_trace_index_t* t, int trc)
{
    DO_SQLITE(sqlite3_bind_int(t->stmt_update_stat, 1, trc));
    DO_INSERT(sqlite3_step(t->stmt_update_stat));
    DO_SQLITE(sqlite3_reset(t->stmt_update_stat));
}

static void segyti_add_trace(segyti_trace_index_t* t, int trc,
                                     double sx, double sy, double gx, double gy)
{
    double off = sqrt((sx-gx)*(sx-gx) + (sy-gy)*(sy-gy));
    double cmpx = 0.5*(sx+gx);
    double cmpy = 0.5*(sy+gy);
    double a = segy_compute_azimuth(sx, sy, gx, gy);

    DO_SQLITE(sqlite3_bind_int(t->stmt_update_xy, 1, trc));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy, 2, sx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy, 3, sx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy, 4, sy));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy, 5, sy));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy, 6, gx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy, 7, gx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy, 8, gy));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy, 9, gy));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy,10, cmpx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy,11, cmpx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy,12, cmpy));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy,13, cmpy));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy,14, a));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy,15, a));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy,16, off));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_xy,17, off));
    DO_INSERT(sqlite3_step(t->stmt_update_xy));
    DO_SQLITE(sqlite3_reset(t->stmt_update_xy));

    if(t->mm.mino > off) t->mm.mino = off;
    if(t->mm.maxo < off) t->mm.maxo = off;

    if(t->mm.mina > a) t->mm.mina = a;
    if(t->mm.maxa < a) t->mm.maxa = a;

    if(t->mm.minsx > sx) t->mm.minsx = sx;
    if(t->mm.maxsx < sx) t->mm.maxsx = sx;
    if(t->mm.minsy > sy) t->mm.minsy = sy;
    if(t->mm.maxsy < sy) t->mm.maxsy = sy;
    if(t->mm.mingx > gx) t->mm.mingx = gx;
    if(t->mm.maxgx < gx) t->mm.maxgx = gx;
    if(t->mm.mingy > gy) t->mm.mingy = gy;
    if(t->mm.maxgy < gy) t->mm.maxgy = gy;

    if(t->mm.mincx > cmpx) t->mm.mincx = cmpx;
    if(t->mm.maxcx < cmpx) t->mm.maxcx = cmpx;
    if(t->mm.mincy > cmpy) t->mm.mincy = cmpy;
    if(t->mm.maxcy < cmpy) t->mm.maxcy = cmpy;
}

static void segyti_update_minmax(segyti_trace_index_t* t)
{
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax,  1, t->mm.mino));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax,  2, t->mm.maxo));

    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax,  3, t->mm.mincx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax,  4, t->mm.maxcx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax,  5, t->mm.mincy));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax,  6, t->mm.maxcy));

    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax,  7, t->mm.minsx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax,  8, t->mm.maxsx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax,  9, t->mm.minsy));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax, 10, t->mm.maxsy));

    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax, 11, t->mm.mingx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax, 12, t->mm.maxgx));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax, 13, t->mm.mingy));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax, 14, t->mm.maxgy));

    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax, 15, t->mm.mina));
    DO_SQLITE(sqlite3_bind_double(t->stmt_update_minmax, 16, t->mm.maxa));

    DO_INSERT(sqlite3_step(t->stmt_update_minmax));
    DO_SQLITE(sqlite3_reset(t->stmt_update_minmax));
}

static void segy_index_create_offset_dist_table(segyti_trace_index_t* t)
{
    if(sqlite3_table_exists(t->db, "offset_dist")) return;
    segyti_dosql(t, "BEGIN TRANSACTION");
    segyti_dosql(t, "CREATE TABLE offset_dist (count INTEGER, o1 REAL)");

    if(sqlite3_table_exists(t->db, "xyoa")) {
        segyti_dosql(t, "INSERT INTO offset_dist SELECT count(*), o1 FROM xyoa GROUP BY(o1)");
    } else {
        segyti_dosql(t, "INSERT INTO offset_dist SELECT count(*), o1 FROM xy GROUP BY(o1)");
    }

    segyti_dosql(t, "COMMIT TRANSACTION");
}

static void segyti_initialize(segyti_trace_index_t* t)
{
    if(t->is_coordinate_file) {
        int64_t size;
        t->sep = sep_open( t->file, SEP_READ, 0 );
        if(t->sep->headers->n[0] != 4) {
            ERROR(("n1=%d in %s: 4 is required for n1 of a coordinates file",
                   t->sep->headers->n[0], t->file));
        }
        size = sep_get_total_size(t->sep);
        t->ntraces = (int)(size/4);
    } else {
        int ns, dt;
        t->segy = segy_get_segy_info( t->file, &ns, &dt, &t->ntraces, t->parset );
    }

    t->nbuf = 655360; /*=10Mb of coordinates: t->nbuf x 4 sx,sy,gx,gy x 4 sieof(float) */
    if(t->nbuf > t->ntraces) t->nbuf = t->ntraces;
    t->coords = alloc1float(t->nbuf*4);

    t->first_idx_in_db = 0;
    t->mm.mino = t->mm.mina = t->mm.minsx = t->mm.mingx = t->mm.mincx =
        t->mm.minsy = t->mm.mingy = t->mm.mincy = 1.0E+30;
    t->mm.maxo = t->mm.maxa = t->mm.maxsx = t->mm.maxgx = t->mm.maxcx =
        t->mm.maxsy = t->mm.maxgy = t->mm.maxcy = -t->mm.mino;

    t->timer_reading = create_timer("Reading coordinates");
    t->timer_updating = create_timer("Updating index");
    t->timer_optimization = create_timer("Optimizing index");

    if(se_io_exists(t->dbfile)) {
        se_io_remove(t->dbfile);
    }

    DO_SQLITE(sqlite3_open_v2(t->dbfile, &t->db,
                              SQLITE_OPEN_READWRITE|SQLITE_OPEN_CREATE,
                              NULL));

    /* segyti_dosql(t, "PRAGMA page_size = 16384"); */
    /* segyti_dosql(t, "PRAGMA cache_size = 4000"); */
    /* PRAGMA journal_mode = DELETE | TRUNCATE | PERSIST | MEMORY | WAL | OFF */
    /* segyti_dosql(t, "PRAGMA journal_mode = OFF");*/
    /* segyti_dosql(t, "PRAGMA read_uncommitted = TRUE"); */

    segyti_dosql(t, "BEGIN TRANSACTION");


    segyti_dosql(t, "CREATE TABLE info (name TEXT PRIMARY KEY, val)");

    segyti_dosql(t,
                 "CREATE TABLE stat (last_trace_recorded INTEGER)");

    segyti_dosql(t,
                 "CREATE TABLE minmax "
                 "(mino REAL, maxo REAL, "
                 " mincx REAL, maxcx REAL, mincy REAL, maxcy REAL,"
                 " minsx REAL, maxsx REAL, minsy REAL, maxsy REAL,"
                 " mingx REAL, maxgx REAL, mingy REAL, maxgy REAL,"
                 " mina REAL, maxa REAL)");

    segyti_dosql(t,
                 "INSERT INTO info VALUES ('version', %d)", SEGYTI_SCHEMA_VERSION);

    segyti_dosql(t,
                 "CREATE VIRTUAL TABLE xy USING rtree"
                 "(id INTEGER, "
                    "sx1 REAL, sx2 REAL, sy1 REAL, sy2 REAL, "
                    "gx1 REAL, gx2 REAL, gy1 REAL, gy2 REAL, "
                    "cx1 REAL, cx2 REAL, cy1 REAL, cy2 REAL, "
                    "a1 REAL, a2 REAL, o1 REAL, o2 REAL)");

    segyti_dosql(t,
                 "INSERT INTO info VALUES ('coordinates_input', %Q)", t->file);
    segyti_dosql(t,
                 "INSERT INTO info VALUES ('first_trace', %d)", t->first_idx_in_db);
    segyti_dosql(t,
                 "INSERT INTO info VALUES ('ox_survey', %f)", t->sd.ox_survey);
    segyti_dosql(t,
                 "INSERT INTO info VALUES ('oy_survey', %f)", t->sd.oy_survey);
    segyti_dosql(t,
                 "INSERT INTO info VALUES ('survey_azimuth', %f)", t->sd.survey_azimuth);
    segyti_dosql(t,
                 "INSERT INTO info VALUES ('left_handed', %d)", t->sd.left_handed);
    segyti_dosql(t,
                 "INSERT INTO info VALUES ('inline_azimuth', %d)", t->sd.inline_azimuth);


    DO_SQLITE(sqlite3_prepare_v2
              ( t->db,
                "INSERT OR REPLACE INTO stat (ROWID, last_trace_recorded) VALUES (1,?)",
                -1, &t->stmt_update_stat, NULL));

    DO_SQLITE(sqlite3_prepare_v2
              ( t->db,
                "INSERT OR REPLACE INTO minmax (ROWID,"
                " mino, maxo,"
                " mincx, maxcx, mincy, maxcy,"
                " minsx, maxsx, minsy, maxsy,"
                " mingx, maxgx, mingy, maxgy, mina, maxa) "
                "VALUES (1,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                -1, &t->stmt_update_minmax, NULL));

    DO_SQLITE(sqlite3_prepare_v2
              ( t->db,
                "INSERT INTO xy VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                -1, &t->stmt_update_xy, NULL));

    segyti_update_stat(t, -1);
    segyti_update_minmax(t);

    segyti_dosql(t, "COMMIT TRANSACTION");
}

static void segyti_destroy(segyti_trace_index_t* t)
{
    DO_SQLITE(sqlite3_finalize(t->stmt_update_stat));
    DO_SQLITE(sqlite3_finalize(t->stmt_update_minmax));
    DO_SQLITE(sqlite3_finalize(t->stmt_update_xy));
    sqlite3_close(t->db);

    if(t->is_coordinate_file) {
        sep_close(t->sep);
    } else {
        segy_destroy_segy_reader(t->segy);
    }

    if(t->coords) {
        ASSERT(t->nbuf > 0);
        free(t->coords);
    }

    destroy_timer(t->timer_optimization);
    destroy_timer(t->timer_updating);
    destroy_timer(t->timer_reading);

    free(t->file);
    if(t->parset) free(t->parset);
    free(t->dbfile);

    if(t->offsets) free(t->offsets);

    free(t);
}

static void segyti_read_coordinates(segyti_trace_index_t* t, int it, int nread)
{
    if(t->is_coordinate_file) {
        se_fsio_seek(t->sep->data->io, (off_t)(it + t->first_idx_in_db)*16);
        se_fsio_read_float(t->sep->data->io, t->coords, 4*nread);
    } else {
        segy_get_segy_traces( t->segy, it + t->first_idx_in_db, nread, NULL, NULL, t->coords, 1);
    }
}

static void segyti_process_coordinates(segyti_trace_index_t* t, int it, int nread)
{
    int i;

    ASSERT(it + nread <= t->ntraces);

    segyti_dosql(t, "BEGIN TRANSACTION");

    for(i = 0; i < nread; ++i) {
        double sx = t->coords[4*i+0];
        double sy = t->coords[4*i+1];
        double gx = t->coords[4*i+2];
        double gy = t->coords[4*i+3];
        to_local_cs( &t->sd, sx, sy, &sx, &sy );
        to_local_cs( &t->sd, gx, gy, &gx, &gy );
        segyti_add_trace(t, it+i+1, sx, sy, gx, gy);
    }

    segyti_update_minmax(t);
    segyti_update_stat(t, it+nread);

    segyti_dosql(t, "COMMIT TRANSACTION");
}

static void segyti_report_timers(segyti_trace_index_t* t, int n)
{
    timer* timers[3] = {t->timer_reading, t->timer_updating, t->timer_optimization};
    report_timers(timers, n);
}

void segy_create_trace_index(const char* file,
                                int is_coordinate_file,
                                const char* outfile,
                                const survey_description_t* sd,
                                const char* parset)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)malloc(sizeof(*t));
    memset(t, 0, sizeof(*t));
    int it, nread;
    time_t last_reported_time = time(NULL);
    int have_optimizations = 0;

    t->file = strdup(file);
    t->parset = parset?strdup(parset):NULL;
    t->dbfile = strdup(outfile);
    t->is_coordinate_file = is_coordinate_file;
    t->sd = *sd;
    t->noffsets = t->offsets_bs = 0;
    t->offsets = NULL;

    segyti_initialize(t);

    INFOV((2, "Creating index file for %d traces with coordinates read from %s", t->ntraces, t->file));

    for(it = 0; it < t->ntraces; it += nread) {
        nread = t->nbuf;
        if(it + nread > t->ntraces) nread = t->ntraces - it;

        start_timer(t->timer_reading);
        segyti_read_coordinates(t, it, nread);
        stop_timer(t->timer_reading);

        INFOV((10, "Indexing %d traces starting at %d", nread, it));

        start_timer(t->timer_updating);
        segyti_process_coordinates(t, it, nread);
        stop_timer(t->timer_updating);
        if(se_get_verb_level() >= 2 && time(NULL) - last_reported_time > 10) {
            last_reported_time = time(NULL);
            segyti_report_timers(t, 2);
        }
    }

    if(se_get_defpar_int32("optimize_oa", 0)) {
        INFOV((3, "Optimizing offset and azimuth queries for the index file %s", t->dbfile));
        start_timer(t->timer_optimization);
        segyti_dosql(t, "BEGIN TRANSACTION");
        segyti_dosql(t, "CREATE VIRTUAL TABLE xyoa USING RTREE(id INTEGER, a1 REAL, a2 REAL, o1 REAL, o2 REAL)");
        segyti_dosql(t, "INSERT INTO xyoa (id,a1,a2,o1,o2) SELECT id,a1,a2,o1,o2 FROM xy");
        segyti_dosql(t, "COMMIT TRANSACTION");
        segy_index_create_offset_dist_table(t);
        stop_timer(t->timer_optimization);
        have_optimizations = 1;
    }

    if(se_get_defpar_int32("optimize_src", 0)) {
        INFOV((3, "Optimizing source queries for the index file %s", t->dbfile));
        start_timer(t->timer_optimization);
        segyti_dosql(t, "BEGIN TRANSACTION");
        segyti_dosql(t, "CREATE VIRTUAL TABLE xys USING RTREE(id INTEGER PRIMARY KEY AUTOINCREMENT, "
                        "sx1 REAL, sx2 REAL, sy1 REAL, sy2 REAL)");
        segyti_dosql(t, "INSERT INTO xys (sx1,sx2,sy1,sy2) SELECT sx1,sx1,sy1,sy1 FROM xy GROUP BY sx1,sy1");
        segyti_dosql(t, "COMMIT TRANSACTION");
        stop_timer(t->timer_optimization);
        have_optimizations = 1;
    }

    if(se_get_defpar_int32("optimize_rec", 0)) {
        INFOV((3, "Optimizing receiver queries for the index file %s", t->dbfile));
        start_timer(t->timer_optimization);
        segyti_dosql(t, "BEGIN TRANSACTION");
        segyti_dosql(t, "CREATE VIRTUAL TABLE xyg USING RTREE(id INTEGER PRIMARY KEY AUTOINCREMENT, "
                        "gx1 REAL, gx2 REAL, gy1 REAL, gy2 REAL)");
        segyti_dosql(t, "INSERT INTO xyg (gx1,gx2,gy1,gy2) SELECT gx1,gx1,gy1,gy1 FROM xy GROUP BY gx1,gy1");
        segyti_dosql(t, "COMMIT TRANSACTION");
        stop_timer(t->timer_optimization);
        have_optimizations = 1;
    }

    if(se_get_defpar_int32("optimize_cmp", 0)) {
        INFOV((3, "Optimizing receiver queries for the index file %s", t->dbfile));
        start_timer(t->timer_optimization);
        segyti_dosql(t, "BEGIN TRANSACTION");
        segyti_dosql(t, "CREATE VIRTUAL TABLE xyc USING RTREE(id INTEGER PRIMARY KEY AUTOINCREMENT, "
                        "cx1 REAL, cx2 REAL, cy1 REAL, cy2 REAL)");
        segyti_dosql(t, "INSERT INTO xyc (cx1,cx2,cy1,cy2) SELECT cx1,cx1,cy1,cy1 FROM xy GROUP BY cx1,cy1");
        segyti_dosql(t, "COMMIT TRANSACTION");
        stop_timer(t->timer_optimization);
        have_optimizations = 1;
    }

    INFOV((2, "Finished creating index file %s for %s", t->dbfile, t->file));
    if(se_get_verb_level() >= 2) {
        segyti_report_timers(t, have_optimizations?3:2);
    }

    segyti_destroy(t);
}


static void segyidx_read_minmax(segyti_trace_index_t* t)
{
    sqlite3_stmt* stmt;
    int col;

    DO_SQLITE(sqlite3_prepare_v2
              ( t->db, "SELECT "
                " mino, maxo,"
                " mincx, maxcx, mincy, maxcy,"
                " minsx, maxsx, minsy, maxsy,"
                " mingx, maxgx, mingy, maxgy, mina, maxa FROM minmax", -1, &stmt, NULL));

    if( sqlite3_step(stmt) != SQLITE_ROW ) {
        ERROR(("Trace index file %s is corrupted: cannot find min/max info in the index file",
               t->dbfile));
    }

    col = -1;
    t->mm.mino = sqlite3_column_double(stmt, ++col);
    t->mm.maxo = sqlite3_column_double(stmt, ++col);

    t->mm.mincx = sqlite3_column_double(stmt, ++col);
    t->mm.maxcx = sqlite3_column_double(stmt, ++col);
    t->mm.mincy = sqlite3_column_double(stmt, ++col);
    t->mm.maxcy = sqlite3_column_double(stmt, ++col);

    t->mm.minsx = sqlite3_column_double(stmt, ++col);
    t->mm.maxsx = sqlite3_column_double(stmt, ++col);
    t->mm.minsy = sqlite3_column_double(stmt, ++col);
    t->mm.maxsy = sqlite3_column_double(stmt, ++col);

    t->mm.mingx = sqlite3_column_double(stmt, ++col);
    t->mm.maxgx = sqlite3_column_double(stmt, ++col);
    t->mm.mingy = sqlite3_column_double(stmt, ++col);
    t->mm.maxgy = sqlite3_column_double(stmt, ++col);

    t->mm.mina = sqlite3_column_double(stmt, ++col);
    t->mm.maxa = sqlite3_column_double(stmt, ++col);

    DO_SQLITE(sqlite3_finalize(stmt));

    INFOV((10,"From index file %s: min/max values are:",t->dbfile));
    INFOV((10,"\toffset  = %g/%g", t->mm.mino, t->mm.maxo));
    INFOV((10,"\tazimuth = %g/%g", t->mm.mina, t->mm.maxa));
    INFOV((10,"\tcmp x,y = %g/%g, %g/%g", t->mm.mincx, t->mm.maxcx, t->mm.mincy, t->mm.maxcy));
    INFOV((10,"\tsrc x,y = %g/%g, %g/%g", t->mm.minsx, t->mm.maxsx, t->mm.minsy, t->mm.maxsy));
    INFOV((10,"\trec x,y = %g/%g, %g/%g", t->mm.mingx, t->mm.maxgx, t->mm.mingy, t->mm.maxgy));
}

static int segyidx_read_info_int(segyti_trace_index_t* t, sqlite3_stmt* stmt, const char* name)
{
    int ret;

    DO_SQLITE(sqlite3_bind_text(stmt, 1, name, -1, SQLITE_STATIC));

    if( sqlite3_step(stmt) != SQLITE_ROW ) {
        ERROR(("Trace index file %s is corrupted: cannot find %s info in the index file",
               t->dbfile, name));
    }

    ret = sqlite3_column_int(stmt, 0);

    DO_SQLITE(sqlite3_reset(stmt));

    return ret;
}

static int segyidx_read_definfo_int(segyti_trace_index_t* t, sqlite3_stmt* stmt, const char* name, int def)
{
    int ret;

    DO_SQLITE(sqlite3_bind_text(stmt, 1, name, -1, SQLITE_STATIC));

    if( sqlite3_step(stmt) == SQLITE_ROW ) {
        ret = sqlite3_column_int(stmt, 0);
    } else {
        ret = def;
    }


    DO_SQLITE(sqlite3_reset(stmt));

    return ret;
}

static double segyidx_read_info_double(segyti_trace_index_t* t, sqlite3_stmt* stmt, const char* name)
{
    double ret;

    DO_SQLITE(sqlite3_bind_text(stmt, 1, name, -1, SQLITE_STATIC));

    if( sqlite3_step(stmt) != SQLITE_ROW ) {
        ERROR(("Trace index file %s is corrupted: cannot find %s info in the index file",
               t->dbfile, name));
    }

    ret = sqlite3_column_double(stmt, 0);

    DO_SQLITE(sqlite3_reset(stmt));

    return ret;
}

static char* segyidx_read_info_str(segyti_trace_index_t* t, sqlite3_stmt* stmt, const char* name)
{
    char* ret;
    const char* text;

    DO_SQLITE(sqlite3_bind_text(stmt, 1, name, -1, SQLITE_STATIC));

    if( sqlite3_step(stmt) != SQLITE_ROW ) {
        ERROR(("Trace index file %s is corrupted: cannot find %s info in the index file",
               t->dbfile, name));
    }
    text = (const char*)sqlite3_column_text(stmt, 0);
    if(text == NULL) {
        text = "";
    }

    ret = strdup(text);

    DO_SQLITE(sqlite3_reset(stmt));

    return ret;
}

static void segyidx_read_survey(segyti_trace_index_t* t, sqlite3_stmt* stmt,
                                        survey_description_t* sd)
{
    init_survey_description(sd);
    sd->ox_survey = segyidx_read_info_double(t, stmt, "ox_survey");
    sd->oy_survey = segyidx_read_info_double(t, stmt, "oy_survey");
    sd->survey_azimuth = segyidx_read_info_double(t, stmt, "survey_azimuth");
    sd->left_handed = segyidx_read_info_int(t, stmt, "left_handed");
    sd->inline_azimuth = segyidx_read_definfo_int(t, stmt, "inline_azimuth", 0);
    survey_description_update_cached_values(&t->sd);
}

static void segyidx_read_info(segyti_trace_index_t* t)
{
    sqlite3_stmt* stmt;
    int version;

    DO_SQLITE(sqlite3_prepare_v2
              ( t->db, "SELECT val FROM info WHERE name=?", -1, &stmt, NULL));

    version = segyidx_read_info_int(t, stmt, "version");
    INFOV((10, "From index file %s: encoding version is %d", t->dbfile, version));
    if(version != SEGYTI_SCHEMA_VERSION) {
        ERROR(("The encoding version %d used for index file %s is different than the one "
               "expected for this version of the program (%d)",
               version, t->dbfile, SEGYTI_SCHEMA_VERSION));
    }

    segyidx_read_survey(t, stmt, &t->sd);
    INFOV((10, "From index file %s: the coordinates transform parameters are:",
           t->dbfile));
    INFOV((10, "\t     ox_survey=%g", t->sd.ox_survey));
    INFOV((10, "\t     oy_survey=%g", t->sd.oy_survey));
    INFOV((10, "\tsurvey_azimuth=%g", t->sd.survey_azimuth));
    INFOV((10, "\t   left_handed=%d", t->sd.left_handed));
    INFOV((10, "\tinline_azimuth=%d", t->sd.inline_azimuth));

    t->first_idx_in_db = segyidx_read_info_int(t, stmt, "first_trace");
    INFOV((10, "From index file %s: first trace index is %d",
           t->dbfile, t->first_idx_in_db));

    t->file = segyidx_read_info_str(t, stmt, "coordinates_input");
    INFOV((10, "From index file %s: original source for coordinates was %s",
           t->dbfile, t->file));

    if(!se_io_exists(t->file)) {
        WARN(("The orignal file %s used to create index file %s doesn't exists.",
              t->file, t->dbfile));
        WARN(("The index may need to be updated or moved to a new location."));
    } else if( se_io_mtime(t->file) > se_io_mtime(t->dbfile) ) {
        WARN(("Index file %s appears to be older than original file %s.",
              t->dbfile, t->file));
        WARN(("The index may need to be updated."));
    }

    DO_SQLITE(sqlite3_finalize(stmt));
}

static void segyidx_check_same_system(segyti_trace_index_t* t,
                                              const survey_description_t* sd)
{
    if( !dequal(t->sd.ox_survey, sd->ox_survey) ||
        !dequal(t->sd.oy_survey, sd->oy_survey) ||
        !dequal(t->sd.survey_azimuth, sd->survey_azimuth) ||
        t->sd.left_handed != sd->left_handed ||
        t->sd.inline_azimuth != sd->inline_azimuth) {
        WARN(("This version doesn't support using the index file with a different "
              "coordinate system than the one used when was generated. "
              "Please re-creeate the index."));
        WARN(("From index file %s: the coordinates transform parameters are:",
              t->dbfile));
        WARN(("\t     ox_survey=%g", t->sd.ox_survey));
        WARN(("\t     oy_survey=%g", t->sd.oy_survey));
        WARN(("\tsurvey_azimuth=%g", t->sd.survey_azimuth));
        WARN(("\t   left_handed=%d", t->sd.left_handed));
        WARN(("\tinline_azimuth=%d", t->sd.inline_azimuth));
        WARN(("Current coordinate transform parameters are:"));
        WARN(("\t     ox_survey=%g", sd->ox_survey));
        WARN(("\t     oy_survey=%g", sd->oy_survey));
        WARN(("\tsurvey_azimuth=%g", sd->survey_azimuth));
        WARN(("\t   left_handed=%d", sd->left_handed));
        WARN(("\tinline_azimuth=%d", sd->inline_azimuth));
        ERROR(("The program will stop now"));
    }
}

inline void segyidx_get_coord_index(sqlite3_context *context, int argc, sqlite3_value **argv, int min)
{
    double sx, sy, gx, gy;
    int res, is, ig;
    segyti_trace_index_t* t;

    ASSERT(argc == 4);
    (void)argc;

    t = (segyti_trace_index_t*)sqlite3_user_data(context);

    sx = sqlite3_value_double(argv[0]);
    sy = sqlite3_value_double(argv[1]);
    gx = sqlite3_value_double(argv[2]);
    gy = sqlite3_value_double(argv[3]);

    is = t->cache_index(t->p_app, sx, sy);
    ig = t->cache_index(t->p_app, gx, gy);
    if(is < ig) res = min?is:ig;
    else        res = min?ig:is;

    sqlite3_result_int(context, res);
}

static void segyidx_min_coord_index(sqlite3_context *context, int argc, sqlite3_value **argv)
{
    segyidx_get_coord_index(context, argc, argv, 1);
}

static void segyidx_max_coord_index(sqlite3_context *context, int argc, sqlite3_value **argv)
{
    segyidx_get_coord_index(context, argc, argv, 0);
}

static void segyidx_filter_trace(sqlite3_context *context, int argc, sqlite3_value **argv)
{
    double sx, sy, gx, gy;
    int res;
    segyti_trace_index_t* t;

    ASSERT(argc == 4);
    (void)argc;

    t = (segyti_trace_index_t*)sqlite3_user_data(context);

    sx = sqlite3_value_double(argv[0]);
    sy = sqlite3_value_double(argv[1]);
    gx = sqlite3_value_double(argv[2]);
    gy = sqlite3_value_double(argv[3]);

    if( t->trace_included )
        res = t->trace_included(t->p_app_check, sx, sy, gx, gy);
    else
        res = 1;

    sqlite3_result_int(context, res);
}

static void segytidx_open(segyti_trace_index_t* t)
{
    char* sql = NULL;
    char* trace_included_str = NULL;
    char* trace_sort_str = NULL;

    if(!se_io_readable_file(t->dbfile)) {
        ERROR(("The index file %s is not readable", t->dbfile));
    }
    DO_SQLITE(sqlite3_open_v2(t->dbfile, &t->db,
                              SQLITE_OPEN_READONLY, NULL));

    segyidx_read_info(t);
    segyidx_read_minmax(t);

    DO_SQLITE(sqlite3_create_function
              (t->db,
               "min_coord_index",
               4, SQLITE_ANY, t,
               segyidx_min_coord_index, NULL, NULL));

    DO_SQLITE(sqlite3_create_function
              (t->db,
               "max_coord_index",
               4, SQLITE_ANY, t,
               segyidx_max_coord_index, NULL, NULL));

    DO_SQLITE(sqlite3_create_function
              (t->db,
               "filter_trace",
               4, SQLITE_ANY, t,
               segyidx_filter_trace, NULL, NULL));

    t->query_idx_minof = 1;
    t->query_idx_maxof = 2;
    t->query_idx_mincx = 3;
    t->query_idx_maxcx = 4;
    t->query_idx_mincy = 5;
    t->query_idx_maxcy = 6;

    if( t->trace_included ) {
        trace_included_str = strdup("filter_trace(sx1,sy1,gx1,gy1) AND ");
    } else {
        trace_included_str = strdup("");
    }

    switch(t->sort_type) {
    case 1:
        INFOV((2, "SEGY Index: Sorting input traces for better cache hit"));
        if(t->cache_index) {
            trace_sort_str = strdup("ORDER by "
                "min_coord_index(sx1,sy1,gx1,gy1),"
                "max_coord_index(sx1,sy1,gx1,gy1)");
        } else {
            trace_sort_str = strdup("ORDER by "
                "MIN(sx1,gx1),MAX(sx1,gx1),MIN(sy1,gy1),MAX(sy1,gy1)");
        }

        break;
    case 2:
        INFOV((2, "SEGY Index: Sorting input traces for better data I/O"));
        trace_sort_str = strdup("ORDER by id");
        break;

    default:
        INFOV((2, "SEGY Index: Not using any special sorting"));
        trace_sort_str = strdup("");
        break;
    }

    sql = asprintf("SELECT id FROM xy WHERE %s"
                       "o1 >= ? AND o1 < ? AND "
                       "cx1 >= ? AND cx1 < ? AND "
                       "cy1 >= ? AND cy1 < ? "
                       "%s", trace_included_str, trace_sort_str);

    DO_SQLITE(sqlite3_prepare_v2
              ( t->db, sql, -1, &t->stmt_get_traces, NULL ));

    t->has_crt_data = 0;
    t->get_traces_binded = 0;
    free(sql);
}

static char* sqlite_db_error(const char* file, sqlite3* db)
{
    char* msg = strdup(sqlite3_errmsg(db));
    char* error = asprintf("Error working with index file %s: %s", file, msg);
    free1char(msg);
    return error;
}

int segy_check_version_trace_index(const char* file, char** error)
{
    int rc, version, success = 0;
    sqlite3* db = NULL;
    sqlite3_stmt* stmt = NULL;

    *error = NULL;

    /* check if the file is readable */
    if(! se_io_readable_file(file)) {
        *error = asprintf("The index file %s is not readable", file);
        return 1;
    }

    while (1) {
        /* open the sqlite database */
        rc = sqlite3_open_v2(file, &db, SQLITE_OPEN_READONLY, NULL);
        if (rc != SQLITE_OK) {
            *error = sqlite_db_error(file, db);
            break;
        }

        /* prepare the statement */
        rc = sqlite3_prepare_v2(db, "SELECT val FROM info WHERE name=?", -1, &stmt, NULL);
        if (rc != SQLITE_OK) {
            *error = sqlite_db_error(file, db);
            break;
        }

        /* bind the arguments */
        rc = sqlite3_bind_text(stmt, 1, "version", -1, SQLITE_STATIC);
        if (rc != SQLITE_OK) {
            *error = sqlite_db_error(file, db);
            break;
        }

        /* execute the statement */
        if (sqlite3_step(stmt) != SQLITE_ROW) {
            *error = asprintf("Trace index file %s is corrupted: cannot find %s info in the index file",
                                  file, "version");
            break;
        }
        version = sqlite3_column_int(stmt, 0);

        /* reset the statement */
        rc = sqlite3_reset(stmt);
        if (rc != SQLITE_OK) {
            *error = sqlite_db_error(file, db);
            break;
        }

        /* validate file version */
        if (version != SEGYTI_SCHEMA_VERSION) {
            *error = asprintf("The encoding version %d used for index file %s is different "
                                  "than the one expected for this version of the program (%d)",
                                  version, file, SEGYTI_SCHEMA_VERSION);
            break;
        }

        success = 1;
        break;
    }

    /* close the statement */
    if (stmt != NULL) {
        rc = sqlite3_finalize(stmt);
        if (success && rc != SQLITE_OK) {
            *error = sqlite_db_error(file, db);
            success = 0;
        }
    }

    /* close the database */
    if (db != NULL) {
        rc = sqlite3_close(db);
        if (success && rc != SQLITE_OK) {
            *error = sqlite_db_error(file, db);
            success = 0;
        }
    }

    return ! success;
}

segy_trace_index_t segy_open_trace_index(const char* file,
                                               const survey_description_t* sd,
                                               void* p_app_check,
                                               int (*trace_included)(void* p, double sx, double sy,
                                                                     double gx, double gy),
                                               void* p_app,
                                               int (*cache_index)(void* p,
                                                                  double x, double y),
                                               int sort_type)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)malloc(sizeof(*t));
    memset(t, 0, sizeof(*t));

    t->noffsets = t->offsets_bs = 0;
    t->offsets = NULL;

    t->dbfile = strdup(file);

    t->p_app = p_app;
    t->cache_index = cache_index;
    t->p_app_check = p_app_check;
    t->trace_included = trace_included;
    t->sort_type = sort_type;

    segytidx_open(t);
    if(sd != NULL)
        segyidx_check_same_system(t, sd);

    return (segy_trace_index_t)t;
}

void segy_close_trace_index(segy_trace_index_t _t)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;

    if(t->offsets) free(t->offsets);

    if(t->has_crt_data)
        DO_SQLITE(sqlite3_reset(t->stmt_get_traces));

    DO_SQLITE(sqlite3_finalize(t->stmt_get_traces));

    sqlite3_close(t->db);

    free(t->file);
    free(t->dbfile);
#ifdef DEBUG
    memset(t, 0, sizeof(*t));
#endif
    free(t);
}

int segy_trace_index_get_traces(segy_trace_index_t _t,
                                   double minof, double maxof,
                                   double mincx, double maxcx,
                                   double mincy, double maxcy,
                                   int* idx, int n)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    int ret = -1, ntr;

    if(!t->get_traces_binded ||
       !dequal(t->last_minof, minof) ||
       !dequal(t->last_maxof, maxof) ||
       !dequal(t->last_mincx, mincx) ||
       !dequal(t->last_maxcx, maxcx) ||
       !dequal(t->last_mincy, mincy) ||
       !dequal(t->last_maxcy, maxcy) ) {
        if(t->has_crt_data)
            DO_SQLITE(sqlite3_reset(t->stmt_get_traces));

        DO_SQLITE(sqlite3_bind_double(t->stmt_get_traces, t->query_idx_minof, minof));
        DO_SQLITE(sqlite3_bind_double(t->stmt_get_traces, t->query_idx_maxof, maxof));
        DO_SQLITE(sqlite3_bind_double(t->stmt_get_traces, t->query_idx_mincx, mincx));
        DO_SQLITE(sqlite3_bind_double(t->stmt_get_traces, t->query_idx_maxcx, maxcx));
        DO_SQLITE(sqlite3_bind_double(t->stmt_get_traces, t->query_idx_mincy, mincy));
        DO_SQLITE(sqlite3_bind_double(t->stmt_get_traces, t->query_idx_maxcy, maxcy));

        t->last_minof = minof;
        t->last_maxof = maxof;
        t->last_mincx = mincx;
        t->last_maxcx = maxcx;
        t->last_mincy = mincy;
        t->last_maxcy = maxcy;

        t->get_traces_binded = 1;
        t->has_crt_data = 1;
    }

    if( ! t->has_crt_data )
        return 0;

    ntr = 0;
    while( ntr < n &&
           (ret = sqlite3_step(t->stmt_get_traces)) == SQLITE_ROW ) {
        idx[ntr++] = sqlite3_column_int(t->stmt_get_traces, 0) - 1;
    }

    t->has_crt_data = (ret == SQLITE_ROW);

    if(!t->has_crt_data)
        DO_SQLITE(sqlite3_reset(t->stmt_get_traces));

    return ntr;
}


static int segy_trace_index_get_traces_any( segyti_trace_index_t* t,
                                                    double minx, double maxx,
                                                    double miny, double maxy,
                                                    int** idx_,
                                                    const char* coordx, const char* coordy )
{
    sqlite3_stmt* stmt;
    int ntr, bufsize;
    int* idx;
    char* sql = asprintf("SELECT id FROM xy "
                             "WHERE %s >= ? AND %s < ? "
                             "AND %s >= ? AND %s < ? ",
                             coordx, coordx, coordy, coordy);


    DO_SQLITE(sqlite3_prepare_v2( t->db, sql , -1, &stmt, NULL ));
    free(sql);

    DO_SQLITE(sqlite3_bind_double(stmt, 1, minx));
    DO_SQLITE(sqlite3_bind_double(stmt, 2, maxx));
    DO_SQLITE(sqlite3_bind_double(stmt, 3, miny));
    DO_SQLITE(sqlite3_bind_double(stmt, 4, maxy));

    ntr = 0;
    bufsize = 0;
    idx = NULL;
    while( sqlite3_step(stmt) == SQLITE_ROW ) {
        if(ntr+1 > bufsize) {
            bufsize += 2048;
            idx = realloc1int(idx, bufsize);
        }
        idx[ntr++] = sqlite3_column_int(stmt, 0) - 1;
    }
    DO_SQLITE(sqlite3_reset(stmt));
    DO_SQLITE(sqlite3_finalize(stmt));

    *idx_ = idx;
    return ntr;
}

inline void segy_min_max_to_global(const survey_description_t* sd,
                                          double* minx, double* maxx,
                                          double* miny, double* maxy)
{
    double x1,y1,x2,y2,x3,y3,x4,y4;
    to_global_cs( sd, *minx, *miny, &x1, &y1);
    to_global_cs( sd, *maxx, *miny, &x2, &y2);
    to_global_cs( sd, *minx, *maxy, &x3, &y3);
    to_global_cs( sd, *maxx, *maxy, &x4, &y4);
    *minx = min4d(x1, x2, x3, x4);
    *maxx = max4d(x1, x2, x3, x4);
    *miny = min4d(y1, y2, y3, y4);
    *maxy = max4d(y1, y2, y3, y4);
}

void segy_trace_index_getminmax(segy_trace_index_t _t, segy_trace_index_statistics_t* st)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    *st = t->mm;

    segy_min_max_to_global(&(t->sd),
                              &(st->minsx), &(st->maxsx), &(st->minsy), &(st->maxsy));
    segy_min_max_to_global(&(t->sd),
                              &(st->mingx), &(st->maxgx), &(st->mingy), &(st->maxgy));
    segy_min_max_to_global(&(t->sd),
                              &(st->mincx), &(st->maxcx), &(st->mincy), &(st->maxcy));
}

int segy_trace_index_get_traces_src( segy_trace_index_t _t,
                                        double minx, double maxx,
                                        double miny, double maxy,
                                        int** idx_ )
{
    return segy_trace_index_get_traces_any((segyti_trace_index_t*)_t,
                                              minx, maxx, miny, maxy, idx_, "sx1", "sy1");
}
int segy_trace_index_get_traces_rec( segy_trace_index_t _t,
                                        double minx, double maxx,
                                        double miny, double maxy,
                                        int** idx_ )
{
    return segy_trace_index_get_traces_any((segyti_trace_index_t*)_t,
                                              minx, maxx, miny, maxy, idx_, "gx1", "gy1");
}

int segy_trace_index_get_noffsets(segy_trace_index_t _t, double delta)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    if(FISZERO(delta) || FISZERO(t->mm.maxo - t->mm.mino)) {
        return 1;
    } else {
        double n = fabs(t->mm.maxo - t->mm.mino) / delta;
        if(is_integer(n)) {
            return (int)n;
        } else {
            return (int)n + 1;
        }
    }
}

double segy_trace_index_get_delta_offset(segy_trace_index_t _t, int n)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    if(n <= 0) n = 1;
    return (t->mm.maxo - t->mm.mino) / n;
}

segy_trace_index_offset_t*
segy_trace_index_get_offset_distribution(segy_trace_index_t _t, int n)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    segy_trace_index_offset_t* od =(segy_trace_index_offset_t*)malloc(n*sizeof(*od));
    memset(od, 0, n*sizeof(*od));
    sqlite3_stmt* stmt;
    int i;
    double d;

    char* sql;
    if(sqlite3_table_exists(t->db, "xyoa")) {
        sql = strdup("SELECT count(*) FROM xyoa WHERE o1 >= ? AND o1 < ?");
    } else {
        sql = strdup("SELECT count(*) FROM xy WHERE o1 >= ? AND o1 < ?");
    }

    DO_SQLITE(sqlite3_prepare_v2( t->db, sql , -1, &stmt, NULL ));
    free(sql);
    d = segy_trace_index_get_delta_offset(_t, n);

    for(i = 0; i < n; ++i) {
        double o0 = t->mm.mino + i*d;
        double o1 = o0 + d + ((i==n-1)?100:0);
        od[i].o = o0;

        DO_SQLITE(sqlite3_bind_double(stmt, 1, o0));
        DO_SQLITE(sqlite3_bind_double(stmt, 2, o1));
        if(sqlite3_step(stmt) != SQLITE_ROW) {
            ERROR(("Internal error accessing the index - expected an offset count"));
        }
        od[i].n = sqlite3_column_int64(stmt, 0);
        DO_SQLITE(sqlite3_reset(stmt));
    }

    DO_SQLITE(sqlite3_finalize(stmt));

    return od;
}

int segy_trace_index_get_nazimuths(segy_trace_index_t _t, double delta)
{
    /*segyti_trace_index_t* t = (segyti_trace_index_t*)_t;*/
    if(FISZERO(delta)) {
        return 1;
    } else {
        double n = 360/delta;/*fabs(t->mm.maxa - t->mm.mina) / delta;*/
        (void)_t;
        if(is_integer(n)) {
            return (int)n;
        } else {
            return (int)n + 1;
        }
    }
}

double segy_trace_index_get_delta_azimuth(segy_trace_index_t _t, int n)
{
    /*segyti_trace_index_t* t = (segyti_trace_index_t*)_t;*/
    (void)_t;
    if(n <= 0) n = 1;
    return 360.0/n; /*(t->mm.maxa - t->mm.mina) / n;*/
}

segy_trace_index_azimuth_t*
segy_trace_index_get_azimuth_distribution(segy_trace_index_t _t, int n)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    segy_trace_index_azimuth_t* ad = (segy_trace_index_azimuth_t*)malloc(n*sizeof(*ad));
    memset(ad, 0, n*sizeof(*ad));

    sqlite3_stmt* stmt;
    int i;
    double d;

    char* sql;
    if(sqlite3_table_exists(t->db, "xyoa")) {
        sql = strdup("SELECT count(*) FROM xyoa WHERE a1 >= ? AND a1 < ?");
    } else {
        sql = strdup("SELECT count(*) FROM xy WHERE a1 >= ? AND a1 < ?");
    }

    DO_SQLITE(sqlite3_prepare_v2( t->db, sql , -1, &stmt, NULL ));
    free(sql);
    d = segy_trace_index_get_delta_azimuth(_t, n);

    for(i = 0; i < n; ++i) {
        double a0 = /*t->mm.mina*/ + i*d;
        double a1 = a0 + d + ((i==n-1)?100:0);
        ad[i].a = a0;

        DO_SQLITE(sqlite3_bind_double(stmt, 1, a0));
        DO_SQLITE(sqlite3_bind_double(stmt, 2, a1));
        if(sqlite3_step(stmt) != SQLITE_ROW) {
            ERROR(("Internal error accessing the index - expected an azimuth count"));
        }
        ad[i].n = sqlite3_column_int64(stmt, 0);
        DO_SQLITE(sqlite3_reset(stmt));
    }

    DO_SQLITE(sqlite3_finalize(stmt));

    return ad;
}

segy_trace_index_oa_t*
segy_trace_index_get_oa_distribution(segy_trace_index_t _t, int no, int na)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    segy_trace_index_oa_t* oad = (segy_trace_index_oa_t*)malloc(no*sizeof(*oad)*na);
    memset(oad, 0, no*sizeof(*oad)*na);


    sqlite3_stmt* stmt;
    int ia, io;
    double dof, daz;

    char* sql;
    if(sqlite3_table_exists(t->db, "xyoa")) {
        sql = strdup("SELECT count(*) FROM xyoa WHERE a1 >= ? AND a1 < ? AND o1 >= ? AND o1 < ?");
    } else {
        sql = strdup("SELECT count(*) FROM xy WHERE a1 >= ? AND a1 < ? AND o1 >= ? AND o1 < ?");
    }
    DO_SQLITE(sqlite3_prepare_v2( t->db, sql , -1, &stmt, NULL ));
    free(sql);
    dof = segy_trace_index_get_delta_offset(_t, no);
    daz = segy_trace_index_get_delta_azimuth(_t, na);

    for(io = 0; io < no; ++io) {
        double o0 = t->mm.mino + io*dof;
        double o1 = o0 + dof + ((io==no-1)?100:0);
        for(ia = 0; ia < na; ++ia) {
            double a0 = /*t->mm.mina*/ + ia*daz;
            double a1 = a0 + daz + ((ia==na-1)?100:0);
            int idx = io*na + ia;
            oad[idx].o = o0;
            oad[idx].a = a0;

            DO_SQLITE(sqlite3_bind_double(stmt, 1, a0));
            DO_SQLITE(sqlite3_bind_double(stmt, 2, a1));
            DO_SQLITE(sqlite3_bind_double(stmt, 3, o0));
            DO_SQLITE(sqlite3_bind_double(stmt, 4, o1));
            if(sqlite3_step(stmt) != SQLITE_ROW) {
                ERROR(("Internal error accessing the index - expected an offset-azimuth count"));
            }
            oad[idx].n = sqlite3_column_int64(stmt, 0);
            DO_SQLITE(sqlite3_reset(stmt));
        }
    }

    DO_SQLITE(sqlite3_finalize(stmt));

    return oad;
}


typedef struct segyti_trace_index_query_s {
    sqlite3_stmt* stmt;
    int have_idx;
    int have_two_points;

    char* sql;
    timer* timer_segy;
    int64_t count;

    int has_stepped;
    segy_trace_index_query_result_t prev;

    const segyti_trace_index_t* t;
} segyti_trace_index_query_t;

static segy_trace_index_query_t
segy_trace_index_query_map(segy_trace_index_t _t,
                              int all,
                              double minx, double maxx,
                              double miny, double maxy,
                              const char* strx, const char* stry,
                              const char* table)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    segyti_trace_index_query_t* tq = ( segyti_trace_index_query_t*)malloc(sizeof(*tq));
    memset(tq, 0, sizeof(*tq));
    tq->t = t;

    {
        double x1,y1,x2,y2,x3,y3,x4,y4;
        to_local_cs( &(t->sd), minx, miny, &x1, &y1);
        to_local_cs( &(t->sd), maxx, miny, &x2, &y2);
        to_local_cs( &(t->sd), minx, maxy, &x3, &y3);
        to_local_cs( &(t->sd), maxx, maxy, &x4, &y4);
        minx = min4d(x1, x2, x3, x4);
        maxx = max4d(x1, x2, x3, x4);
        miny = min4d(y1, y2, y3, y4);
        maxy = max4d(y1, y2, y3, y4);
    }

    if(sqlite3_table_exists(t->db, table)) {
        if(all) {
            tq->sql = asprintf("SELECT %s, %s FROM %s", strx, stry, table);
        } else {
            tq->sql = asprintf("SELECT %s, %s FROM %s WHERE %s >= ? AND %s < ? AND %s >= ? AND %s < ?",
                                   strx, stry, table, strx, strx, stry, stry);
        }
    } else {
        if(all) {
            tq->sql = asprintf("SELECT DISTINCT %s, %s FROM xy", strx, stry);
        } else {
            tq->sql = asprintf("SELECT DISTINCT %s, %s FROM xy WHERE %s >= ? AND %s < ? AND %s >= ? AND %s < ?",
                                   strx, stry, strx, strx, stry, stry);
        }
    }

    if(se_get_verb_level() >= 100) {
        tq->timer_segy = create_timer("query");
        INFOV((500, "Initating map query [%s]", tq->sql));
    }

    DO_SQLITE(sqlite3_prepare_v2( t->db, tq->sql , -1, &tq->stmt, NULL ));

    if(!all) {
        DO_SQLITE(sqlite3_bind_double(tq->stmt, 1, minx));
        DO_SQLITE(sqlite3_bind_double(tq->stmt, 2, maxx));
        DO_SQLITE(sqlite3_bind_double(tq->stmt, 3, miny));
        DO_SQLITE(sqlite3_bind_double(tq->stmt, 4, maxy));
    }

    tq->have_idx = 0;
    tq->have_two_points = 0;
    tq->has_stepped = 0;

    return (segy_trace_index_query_t)tq;
}

void segy_trace_index_destroy_query(segy_trace_index_query_t _q)
{
    segyti_trace_index_query_t* tq = (segyti_trace_index_query_t*)_q;
    sqlite3_finalize(tq->stmt);
    if(se_get_verb_level() >= 100) {
        INFOV((500, "Destroying query [%s]", tq->sql));
        if(tq->sql) free(tq->sql);
        if(tq->timer_segy) destroy_timer(tq->timer_segy);
    }

    tq->stmt = NULL;
    tq->t = NULL;
    free(tq);
}

static int segy_trace_index_same_result(const segyti_trace_index_query_t* tq,
                                                const segy_trace_index_query_result_t* res)
{
    if(tq->has_stepped) {
        if(tq->have_idx && res->id != tq->prev.id) return 0;
        if(FNOTEQUAL(res->x, tq->prev.x) || FNOTEQUAL(res->y, tq->prev.y)) return 0;
        if(tq->have_two_points) {
            if(FNOTEQUAL(res->x1, tq->prev.x1) || FNOTEQUAL(res->y1, tq->prev.y1)) return 0;
        }
        return 1;
    }

    return 0;
}

int segy_trace_index_query_next(segy_trace_index_query_t _q,
                                   segy_trace_index_query_result_t* res)
{
    segyti_trace_index_query_t* tq = (segyti_trace_index_query_t*)_q;

    if(!tq->has_stepped)
        tq->count = 0;

    if(tq->timer_segy)
        start_timer(tq->timer_segy);

    for(;;) {
        int col;

        if(sqlite3_step(tq->stmt) != SQLITE_ROW) {
            sqlite3_reset(tq->stmt);
            if(tq->timer_segy) {
                stop_timer(tq->timer_segy);
                INFOV((100, "Finished query: %Ld items in %gs/%gs; sql:[%s]",
                       (long long int)tq->count,
                       tq->timer_segy->total_atime, tq->timer_segy->total_etime, tq->sql));
            }
            return 0;
        }
        ++tq->count;

        col = 0;
        if(tq->have_idx) {
            res->id = sqlite3_column_int64(tq->stmt, col++);
        }

        res->x = sqlite3_column_double(tq->stmt, col++);
        res->y = sqlite3_column_double(tq->stmt, col++);
        to_global_cs( &(tq->t->sd), res->x, res->y, &(res->x), &(res->y));

        if(tq->have_two_points) {
            res->x1 = sqlite3_column_double(tq->stmt, col++);
            res->y1 = sqlite3_column_double(tq->stmt, col++);
            to_global_cs( &(tq->t->sd), res->x1, res->y1, &(res->x1), &(res->y1));
        }

        if(!segy_trace_index_same_result(tq, res)) {
            tq->has_stepped = 1;
            tq->prev = *res;
            break;
        }

        tq->has_stepped = 1;
    }

    if(tq->timer_segy)
        stop_timer(tq->timer_segy);

    return 1;
}

segy_trace_index_query_t segy_trace_index_query_src(segy_trace_index_t _t,
                                                          int all,
                                                          double minx, double maxx,
                                                          double miny, double maxy)
{
    return segy_trace_index_query_map(_t, all, minx, maxx, miny, maxy, "sx1", "sy1", "xys");
}

segy_trace_index_query_t segy_trace_index_query_rec(segy_trace_index_t _t,
                                                          int all,
                                                          double minx, double maxx,
                                                          double miny, double maxy)
{
    return segy_trace_index_query_map(_t, all, minx, maxx, miny, maxy, "gx1", "gy1", "xyg");
}

segy_trace_index_query_t segy_trace_index_query_cmp(segy_trace_index_t _t,
                                                          int all,
                                                          double minx, double maxx,
                                                          double miny, double maxy)
{
    return segy_trace_index_query_map(_t, all, minx, maxx, miny, maxy, "cx1", "cy1", "xyc");
}


static segy_trace_index_query_t
segy_trace_index_query_gather(segy_trace_index_t _t,
                                 double x, double y, double dx, double dy,
                                 const char* strxy,
                                 const char* strx, const char* stry,
                                 int have_two_points)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    segyti_trace_index_query_t* tq = (segyti_trace_index_query_t*)malloc(sizeof(*tq));
    memset(tq, 0, sizeof(*tq));

    tq->t = t;

    to_local_cs( &(t->sd), x, y, &x, &y);

    tq->sql = asprintf("SELECT id, %s FROM xy WHERE %s>=? AND %s<=? AND %s>=? AND %s<=?",
                       strxy, strx, strx, stry, stry);
    if(se_get_verb_level() >= 100) {
        tq->timer_segy = create_timer("query");
        INFOV((500, "Initating gather query [%s]; x=%f, y=%f, dx=%f, dy=%f",
               tq->sql, x, y, dx, dy));
    }

    DO_SQLITE(sqlite3_prepare_v2( t->db, tq->sql , -1, &tq->stmt, NULL ));

    dx = fabs(0.5*dx);
    dy = fabs(0.5*dy);

    DO_SQLITE(sqlite3_bind_double(tq->stmt, 1, x-dx));
    DO_SQLITE(sqlite3_bind_double(tq->stmt, 2, x+dx));
    DO_SQLITE(sqlite3_bind_double(tq->stmt, 3, y-dy));
    DO_SQLITE(sqlite3_bind_double(tq->stmt, 4, y+dy));

    tq->have_idx = 1;
    tq->have_two_points = have_two_points;
    tq->has_stepped = 0;

    return (segy_trace_index_query_t)tq;
}

segy_trace_index_query_t
segy_trace_index_query_gather_src(segy_trace_index_t _t,
                                     double x, double y, double dx, double dy)
{
    return segy_trace_index_query_gather(_t, x, y, dx, dy, "gx1, gy1", "sx1", "sy1", 0);
}

segy_trace_index_query_t
segy_trace_index_query_gather_rec(segy_trace_index_t _t,
                                     double x, double y, double dx, double dy)
{
    return segy_trace_index_query_gather(_t, x, y, dx, dy, "sx1, sy1", "gx1", "gy1", 0);
}

segy_trace_index_query_t
segy_trace_index_query_gather_cmp(segy_trace_index_t _t,
                                     double x, double y, double dx, double dy)
{
    return segy_trace_index_query_gather(_t, x, y, dx, dy, "sx1, sy1, gx1, gy1", "cx1", "cy1", 1);
}


static int segy_trace_index_get_close_trace(segy_trace_index_t _t,
                                                    double x, double y, double dx, double dy,
                                                    segy_trace_index_query_result_t* res,
                                                    const char* strxy,
                                                    const char* strx, const char* stry,
                                                    const char* table)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    sqlite3_stmt* stmt;
    char* sql;
    int ret;
    double mind;
    timer* timer = NULL;
    int64_t count;

    to_local_cs( &(t->sd), x, y, &x, &y);

    dx = fabs(0.5*dx);
    dy = fabs(0.5*dy);

    if(sqlite3_table_exists(t->db, table)) {
        sql = asprintf("SELECT id, %s FROM %s WHERE %s>=? AND %s<=? AND %s>=? AND %s<=?",
                           strxy, table, strx, strx, stry, stry);
    } else {
        sql = asprintf("SELECT id, %s FROM xy WHERE %s>=? AND %s<=? AND %s>=? AND %s<=?",
                           strxy, strx, strx, stry, stry);
    }

    DO_SQLITE(sqlite3_prepare_v2( t->db, sql, -1, &stmt, NULL ));
    if(se_get_verb_level() >= 100) {
        INFOV((500, "Doing [%s]", sql));
        timer = create_timer("query");
        start_timer(timer);
    }

    DO_SQLITE(sqlite3_bind_double(stmt, 1, x-dx));
    DO_SQLITE(sqlite3_bind_double(stmt, 2, x+dx));
    DO_SQLITE(sqlite3_bind_double(stmt, 3, y-dy));
    DO_SQLITE(sqlite3_bind_double(stmt, 4, y+dy));

    ret = 0;
    mind = 1.0E+20;
    count = 0;
    while(sqlite3_step(stmt) == SQLITE_ROW) {
        double x1 = sqlite3_column_double(stmt, 1);
        double y1 = sqlite3_column_double(stmt, 2);
        double d = (x-x1)*(x-x1) + (y-y1)*(y-y1);
        if(mind > d || !ret) {
            res->id = sqlite3_column_int64(stmt, 0);
            to_global_cs( &(t->sd), x1, y1, &x1, &y1);
            res->x = x1;
            res->y = y1;
            mind = d;
        }
        ret = 1;
        ++count;
    }

    sqlite3_finalize(stmt);
    if(se_get_verb_level() >= 100) {
        stop_timer(timer);
        to_global_cs( &(t->sd), x, y, &x, &y);
        INFOV((100, "Finished close trace query: %Ld items in %gs/%gs; (%g,%g +/- %g,%g); sql:[%s]",
               (long long int)count, timer->total_atime, timer->total_etime,
               x, y, dx, dy,
               sql));
        destroy_timer(timer);
    }

    free(sql);

    return ret;
}

int segy_trace_index_get_close_src(segy_trace_index_t _t,
                                      double x, double y, double dx, double dy,
                                      segy_trace_index_query_result_t* res)
{
    return segy_trace_index_get_close_trace(_t, x, y, dx, dy, res, "sx1, sy1", "sx1", "sy1", "xys");
}

int segy_trace_index_get_close_rec(segy_trace_index_t _t,
                                      double x, double y, double dx, double dy,
                                      segy_trace_index_query_result_t* res)
{
    return segy_trace_index_get_close_trace(_t, x, y, dx, dy, res, "gx1, gy1", "gx1", "gy1", "xyg");
}

int segy_trace_index_get_close_cmp(segy_trace_index_t _t,
                                      double x, double y, double dx, double dy,
                                      segy_trace_index_query_result_t* res)
{
    return segy_trace_index_get_close_trace(_t, x, y, dx, dy, res, "cx1, cy1", "cx1", "cy1", "xyc");
}


char* segy_get_index_file(const char* segy_par_name)
{
    char* tmp;
    char* dir;
    char* ifile = NULL;

    if(segy_par_name == NULL) segy_par_name = "data";

    if(se_have_namedpar(segy_par_name, "index_file", 0)) {
        ifile = se_get_namedpar_str(segy_par_name, "index_file");
        if(se_is_path_absolute(ifile))
            return ifile;
    }

    if( se_have_par(segy_par_name) ) {
        char* dfile;
        se_jcs_set_file_parameter_to_jobdir(segy_par_name);
        dfile = se_get_par_str(segy_par_name);
        dir = se_get_base_path(dfile);
        if(!ifile)
            ifile = asprintf("%s.idx", dfile);
        free(dfile);
    } else {
        dir = se_jcs_get_jobdir_copy();
    }
    if(!ifile)
        ifile = strdup("index_file.idx");

    tmp = se_make_path_absolute( ifile, dir );
    free(dir);
    free(ifile);
    return tmp;
}

char* segy_get_coordinates_file(const char* segy_par_name)
{
    char* tmp;
    char* dir;
    char* ifile = NULL;

    if(segy_par_name == NULL) segy_par_name = "data";

    if(se_have_namedpar(segy_par_name, "coordinates_file", 0)) {
        ifile = se_get_namedpar_str(segy_par_name, "coordinates_file");
        if(se_is_path_absolute(ifile))
            return ifile;
    } else if(se_have_par("coordinates_file")) {
        ifile = se_get_par_str("coordinates_file");
        if(se_is_path_absolute(ifile))
            return ifile;
    }

    if( se_have_par(segy_par_name) ) {
        char* dfile;
        se_jcs_set_file_parameter_to_jobdir(segy_par_name);
        dfile = se_get_par_str(segy_par_name);
        dir = se_get_base_path(dfile);
        if(!ifile)
            ifile = asprintf("%s.coordinates.H", dfile);
        free(dfile);
    } else {
        dir = se_jcs_get_jobdir_copy();
    }
    if(!ifile)
        ifile = strdup("coordinates_file.H");

    tmp = se_make_path_absolute( ifile, dir );
    free(dir);
    free(ifile);
    return tmp;
}

static void segy_index_ensure_offset_distribution_loaded(segyti_trace_index_t* t)
{
    const char* sql;
    sqlite3_stmt* stmt;
    int grow_by = 8*1024;

    if(t->offsets) return;
    t->offsets_bs = 0;
    t->noffsets = 0;

    if(!sqlite3_table_exists(t->db, "offset_dist") && !sqlite3_db_readonly(t->db, t->dbfile)) {
        segy_index_create_offset_dist_table(t);
        sql = "SELECT count, o1 FROM offset_dist ORDER BY o1";
    } else {
        if(sqlite3_table_exists(t->db, "xyoa")) {
            sql = "SELECT count(*), o1 FROM xyoa GROUP BY(o1) ORDER BY o1";
        } else {
            sql = "SELECT count(*), o1 FROM xy GROUP BY(o1) ORDER BY o1";
        }
    }

    DO_SQLITE(sqlite3_prepare_v2( t->db, sql, -1, &stmt, NULL ));
    INFOV((500, "Doing [%s]", sql));
    while(sqlite3_step(stmt) == SQLITE_ROW) {
        if(t->offsets_bs < t->noffsets+1) {
            t->offsets_bs = t->noffsets+1+grow_by;
            t->offsets = (segy_index_offset_bin_t*)realloc(t->offsets, t->offsets_bs*sizeof(*(t->offsets)));
        }
        t->offsets[t->noffsets].cnt = sqlite3_column_int(stmt, 0);
        t->offsets[t->noffsets].off = (float)sqlite3_column_double(stmt, 1);
        t->noffsets += 1;
    }
    sqlite3_finalize(stmt);
}

int segy_index_compute_offset_distribution(segy_trace_index_t _t,
                                              double min, double max, int nbins,
                                              double **values,
                                              int **counts,
                                              int* nvalues)
{
    segyti_trace_index_t* t = (segyti_trace_index_t*)_t;
    int i, i0, i1, ntot;

    segy_index_ensure_offset_distribution_loaded(t);
    if(fabs(min - t->offsets[0].off) < 0.001) min -= fabs(min);
    if(fabs(max - t->offsets[t->noffsets-1].off) < 0.001) max += fabs(max);

    for(i0 = 0; i0 < (int)t->noffsets && t->offsets[i0].off < min; ++i0);
    for(i1 = i0, ntot = 0; i1 < (int)t->noffsets && t->offsets[i1].off <= max; ++i1)
        ntot += t->offsets[i1].cnt;
    INFOV((500, "Found %d traces with offset between %g (%d) and %g (%d); there are %Ld total offset values",
           ntot, min, i0, max, i1, (long long int)t->noffsets));

    *nvalues = nbins + 1;
    *values = alloc1double(*nvalues);
    if(counts != NULL)
        *counts = alloc1int_zero(nbins);

    if(ntot == 0) {
        for(i = 0; i < nbins; ++i) {
            (*values)[i] = min + i*(max-min)/nbins;
        }
        (*values)[nbins] = max;
    } else {
        int ninbin, cnt, crt;
        ninbin = ntot/nbins;
        if(ninbin < 1) ninbin = 1;

        INFOV((500, "Will try to fit about %d traces per bin", ninbin));

        (*values)[0] = t->offsets[i0].off;
        for(crt = 0, cnt = 0, i = i0; i < i1; ++i) {
            cnt += t->offsets[i].cnt;

            if(counts) (*counts)[crt] = cnt;
            if(i < i1-1) (*values)[crt+1] = t->offsets[i+1].off;
            else         (*values)[crt+1] = t->offsets[i].off;

            if(cnt >= ninbin && crt < nbins-1 && i < i1-1) {
                cnt = 0;
                ++crt;
            }
        }
        *nvalues = crt+2;
    }

    return CODE_SUCCESS;
}
