#include "../include/se_par_sep.h"
#include "../include/se_fs_segy.h"
#include "../include/se_fs_sep.h"
#include "../include/se_fs_segy_trchead_eval.h"

typedef struct segy_file_s {
    char* file_name;
    byte reel[EBC_BYTES + BHEAD_BYTES];
    int format, bps, ns, idt;
    int ntraces;
    se_fsio* io_data;
    byte* trace; /* pointer to raw trace, including header header;
                    at some point is cast to ztr_trchead* */
    float* trace_data;
    /* function to convert input data to float */
    segy_fix_data_fnc fix_data;
    segy_ieee2segy_fnc ieee2segy;
} segy_file;

typedef struct extract_traces_s {
    segy_file** data;
    int nfiles;
    segy_file *last_sf_open;

    segy_reader_params_t rp;
} extract_traces;

static extract_traces* init_segy_read(const char* file_name,
                                      const segy_reader_params_t* rp,
                                      const char* parset);
static segy_file* read_one_trace(extract_traces* et, int t, int header_only, int raw);


int segy_get_segy_traces( segy_read_t segy,
                             int first_trace, int ntraces, const int* index_traces,
                             float* data, float* coords, int coords_only)
{
    int nt;
    extract_traces* et = (extract_traces*)segy;

    if(et->nfiles == 0) ERROR(("No files found"));

    if( coords_only == 1 ) INFOV((20, "Reading coordinates only"));

    for(nt = 0; nt < ntraces; nt++) {
        segy_file* sf;
        float* d = data + nt*et->rp.ns;
        int itrace;

        if(index_traces != NULL)
            itrace = index_traces[nt];
        else
            itrace = first_trace+nt;

        sf = read_one_trace(et, itrace, coords_only, 0);

        if(sf == NULL){
            WARN(("Not able to read trace %d from input segy", itrace));
            break;
        }
        if (coords_only != 1) {
            if(sf->ns >= et->rp.ns) {
                memcpy(d, sf->trace_data, et->rp.ns*sizeof(float));
            } else {
                memcpy(d, sf->trace_data, sf->ns*sizeof(float));
                int i;
                for(i = sf->ns; i < et->rp.ns; ++i) {
                    d[i] = 0.0f;
                }
            }

            if(et->rp.sepweights) {
                int i;
                float w;
                se_fsio_seek(et->rp.sepweights->data->io, (off_t)itrace*4);
                se_fsio_read_float(et->rp.sepweights->data->io, &w, 1);
                for(i = 0; i < et->rp.ns; ++i) {
                    d[i] *= w;
                }
            }
        }

        if(coords != NULL)
            segy_compute_trace_coordinates( et, (trchead*)sf->trace,
                                               coords + nt*4 + 0,
                                               coords + nt*4 + 1,
                                               coords + nt*4 + 2,
                                               coords + nt*4 + 3 );
    }

    return nt;
}

int segy_get_segy_traces_and_headers( segy_read_t segy,
                                         int first_trace,
                                         int ntraces, const int* index_traces,
                                         float* data, float* coords,
                                         trchead* hdrs )
{
    int nt;
    extract_traces* et = (extract_traces*)segy;

    if(et->nfiles == 0) ERROR(("No files found"));

    if( data == NULL ) INFOV((10, "Reading headers only"));

    for(nt = 0; nt < ntraces; nt++) {
        segy_file* sf;
        float* d = data + nt*et->rp.ns;
        float* c = coords + nt*4;
        int itrace;

        if(index_traces != NULL)
            itrace = index_traces[nt];
        else
            itrace = first_trace+nt;

        sf = read_one_trace(et, itrace, data==NULL, 1);
        
        if(sf == NULL){
            WARN(("Not able to read trace %d from input segy", itrace));
            break;
        }

        segy_trcheader_order_bytes((trchead*)sf->trace);

        if(hdrs) {
            memcpy(hdrs + nt, sf->trace, THEAD_BYTES);
        }

        if(coords) {
            segy_compute_trace_coordinates( et, (trchead*)sf->trace,
                                               c + 0, c + 1, c + 2, c + 3 );
        }

        if( data ) {
            sf->fix_data((trchead*)sf->trace,
                         sf->trace + THEAD_BYTES, sf->trace_data, sf->ns);

            if(sf->ns >= et->rp.ns) {
                memcpy(d, sf->trace_data, et->rp.ns*sizeof(float));
            } else {
                memcpy(d, sf->trace_data, sf->ns*sizeof(float));
                int i;
                for(i = sf->ns; i < et->rp.ns; ++i) {
                    d[i] = 0.0f;
                }
            }

            if(et->rp.sepweights) {
                int i;
                float w;
                se_fsio_seek(et->rp.sepweights->data->io, (off_t)itrace*4);
                se_fsio_read_float(et->rp.sepweights->data->io, &w, 1);
                for(i = 0; i < et->rp.ns; ++i) {
                    d[i] *= w;
                }
            }
        }
    }

    return nt;
}

int segy_get_raw_segy_traces( segy_read_t segy,
                                 int first_trace, int ntraces,
                                 const int* index_traces,
                                 byte* data)
{
    int nt;
    extract_traces* et = (extract_traces*)segy;

    if(et->nfiles == 0) ERROR(("No files found"));


    for(nt = 0; nt < ntraces; nt++) {
        segy_file* sf;
        byte* d;
        int itrace;

        if(index_traces != NULL)
            itrace = index_traces[nt];
        else
            itrace = first_trace+nt;

        sf = read_one_trace(et, itrace, 0, 1);

        if(sf == NULL){
            WARN(("Not able to read trace %d from input segy",
                  (index_traces)?index_traces[nt]:first_trace+nt ));
            break;
        }

        d = data + (size_t)nt*(THEAD_BYTES + et->rp.ns*sf->bps);

        if(sf->ns >= et->rp.ns) {
            memcpy(d, sf->trace, THEAD_BYTES + et->rp.ns*sf->bps);
        } else {
            memcpy(d, sf->trace, THEAD_BYTES + sf->ns*sf->bps);
            memset(d + THEAD_BYTES + sf->ns*sf->bps, 0, (et->rp.ns-sf->ns)*sf->bps);
        }
        if(sf->ns != et->rp.ns) {
            segy_trcheader_order_bytes((trchead*)d);
            ((trchead*)d)->ns = (unsigned short)et->rp.ns;
            segy_trcheader_order_bytes((trchead*)d);
        }
    }

    return nt;
}


int segy_get_segy_bps( const segy_read_t segy )
{
    const extract_traces* et = (const extract_traces*)segy;
    const segy_file *sf = et->data[0];
    return sf->bps;
}

int segy_get_trace_length( const segy_read_t segy )
{
    const extract_traces* et = (const extract_traces*)segy;
    const segy_file *sf = et->data[0];
    return THEAD_BYTES + et->rp.ns*sf->bps;
}

const byte* segy_get_reel_header(const segy_read_t segy)
{
    const extract_traces* et = (const extract_traces*)segy;
    const segy_file *sf = et->data[0];
    return sf->reel;
}

char* segy_get_ebcdic_text(const segy_read_t segy)
{
    const extract_traces* et = (const extract_traces*)segy;
    const segy_file *sf = et->data[0];
    char* ebcdic = alloc1char(EBC_BYTES);
    memcpy(ebcdic, sf->reel, EBC_BYTES);
    ebcdic2ascii( ebcdic, EBC_BYTES );
    return ebcdic;
}

void segy_get_bhead(const segy_read_t segy, binhead* bhead)
{
    const extract_traces* et = (const extract_traces*)segy;
    const segy_file *sf = et->data[0];

    memcpy(bhead, sf->reel+EBC_BYTES, BHEAD_BYTES);
    segy_binheader_order_bytes(bhead);
}

int segy_get_ns(const segy_read_t segy)
{
    const extract_traces* et = (const extract_traces*)segy;
    return et->rp.ns;
}

int segy_get_ntraces(const segy_read_t segy)
{
    const extract_traces* et = (const extract_traces*)segy;
    int i, ntraces = 0;
    for(i = 0; i < et->nfiles; i++) {
        ntraces += et->data[i]->ntraces;
    }
    return ntraces;
}
int segy_get_dt(const segy_read_t segy)
{
    const extract_traces* et = (const extract_traces*)segy;
    return et->data[0]->idt;
}
int segy_get_format(const segy_read_t segy)
{
    const extract_traces* et = (const extract_traces*)segy;
    return et->data[0]->format;
}
int segy_get_nfiles(const segy_read_t segy)
{
    const extract_traces* et = (const extract_traces*)segy;
    return et->nfiles;
}

segy_read_t segy_create_reader( const char* file_name,
                                      segy_reader_params_t* p)
{
    extract_traces* et;

    et = init_segy_read(file_name, p, NULL);
    return (segy_read_t)et;
}

segy_read_t segy_get_segy_info( const char* file_name,
                                      int* ns, int* dt,
                                      int* ntraces,
                                      const char* parset)
{
    extract_traces* et;

    et = init_segy_read(file_name, NULL, parset);
    if(et->nfiles == 0) ERROR(("No files found"));

    *ns = et->rp.ns;
    *dt = et->data[0]->idt;

    *ntraces = segy_get_ntraces((segy_read_t)et);

    return (segy_read_t)et;
}


static segy_file* read_one_trace( extract_traces* et,
                                  const int t, int header_only, int raw )
{
    int nt, file;

    for( nt = 0, file = 0; file < et->nfiles; file++) {
        segy_file *sf = et->data[file];
        nt += sf->ntraces;
        if(nt > t) {
            int filet = t - (nt - sf->ntraces);
            off_t offset = EBC_BYTES + BHEAD_BYTES +
                (off_t)filet * (THEAD_BYTES + sf->ns*sf->bps);
            INFOV((510, "Found trace %d in %s:%d", t+1,
                   sf->file_name, filet+1));

            if(!sf->io_data) {
                if(et->last_sf_open && !et->rp.keep_open) {
                    if(et->last_sf_open->io_data)
                        se_fsio_close(et->last_sf_open->io_data);
                    et->last_sf_open->io_data = NULL;
                }
                INFOV((500, "Need to open file %s", sf->file_name));
                sf->io_data = se_fsio_init(sf->file_name, "r", NULL);
                et->last_sf_open = sf;
                if(et->rp.use_mmap) {
                    se_fsio_set_mmap(sf->io_data, et->rp.mmap_flags);
                }
            }

            se_fsio_seek(sf->io_data, offset);
            if( header_only == 0 )
                se_fsio_raw_read(sf->io_data, sf->trace,
                                THEAD_BYTES + sf->ns*sf->bps);
            else
                se_fsio_raw_read(sf->io_data, sf->trace, THEAD_BYTES );

            if(!raw) {
                segy_trcheader_order_bytes((trchead*)sf->trace);
                sf->fix_data((trchead*)sf->trace,
                             sf->trace + THEAD_BYTES, sf->trace_data, sf->ns);
            }

            return sf;
        }
    }

    return NULL;
}


void segy_init_reader_par_from_parset(segy_reader_params_t* p, const char* parset)
{
    segy_init_reader_par(p);

    if( ! se_par_initialized() ) return;

    if(parset) p->forced_format = se_get_defnamedpar_int32(parset, "forced_format", p->forced_format, 0);
    else       p->forced_format = se_get_defpar_int32("forced_format", p->forced_format);

    if(parset) p->forced_idt = se_get_defnamedpar_int32(parset, "forced_dt", p->forced_idt, 0);
    else       p->forced_idt = se_get_defpar_int32("forced_dt", p->forced_idt);

    if(parset) p->forced_ns = se_get_defnamedpar_int32(parset, "forced_ns", p->forced_ns, 0);
    else       p->forced_ns = se_get_defpar_int32("forced_ns", p->forced_ns);

    if(parset) p->ns = se_get_defnamedpar_int32(parset, "ns", p->ns, 0);
    else       p->ns = se_get_defpar_int32("ns", p->ns);

    if(parset) p->use_scalco = se_get_defnamedpar_int32(parset, "use_scalco", p->use_scalco, 0);
    else       p->use_scalco = se_get_defpar_int32("use_scalco", p->use_scalco);

    if(parset) p->scalco_factor = se_get_defnamedpar_double(parset, "scalco_factor", p->scalco_factor, 0);
    else       p->scalco_factor = se_get_defpar_double("scalco_factor", p->scalco_factor);

#define DEAL_WITH_EXPRESSION(pr, p_str)                                 \
    if( parset && se_have_namedpar(parset, #pr, 0) ) {                 \
        p->p_str = se_get_namedpar_str(parset, #pr);                   \
    } else if( se_have_par(#pr) ) {                                    \
        p->p_str = se_get_par_str(#pr);                                \
    }

    DEAL_WITH_EXPRESSION(sx, sx_str);
    DEAL_WITH_EXPRESSION(sy, sy_str);
    DEAL_WITH_EXPRESSION(gx, gx_str);
    DEAL_WITH_EXPRESSION(gy, gy_str);

    segy_init_reader_par_expressions(p);

    if( parset && se_have_namedpar(parset, "weights", 0) ) {
        p->weights = se_get_namedpar_str(parset, "weights");
    } else if( se_have_par("weights") ) {
        p->weights = se_get_par_str("weights");
    }

    if(parset) p->keep_open = se_get_defnamedpar_int32(parset, "keep_open", p->keep_open, 0);
    else       p->keep_open = se_get_defpar_int32("keep_open", p->keep_open);

    if(parset) p->use_mmap = se_get_defnamedpar_int32(parset, "mmap", p->use_mmap, 0);
    else       p->use_mmap = se_get_defpar_int32("mmap", p->use_mmap);
    if( parset && se_have_namedpar(parset, "mmap_flags", 0) ) {
        p->mmap_flags = se_get_namedpar_str(parset, "mmap_flags");
    } else if( se_have_par("mmap_flags") ) {
        p->mmap_flags = se_get_par_str("mmap_flags");
    }
}

void segy_init_reader_par_expression_strings(segy_reader_params_t* p,
                                                const char* sx,
                                                const char* sy,
                                                const char* gx,
                                                const char* gy)
{
    if(sx) p->sx = se_strdup(sx);
    if(sy) p->sy = se_strdup(sy);
    if(gx) p->gx = se_strdup(gx);
    if(gy) p->gy = se_strdup(gy);
    segy_init_reader_par_expressions(p);
}

void segy_init_reader_par_expressions(segy_reader_params_t* p)
{
#define INIT_EXPRESSION(pr, p_str)                                      \
    if( p->p_str ) {                                                    \
        p->pr = init_trchead_eval( p->p_str );                      \
        if( p->pr == NULL )                                             \
            WARN(("Problems with %s replacement expression", #pr));     \
    } else {                                                            \
        p->pr = NULL;                                                   \
    }

    INIT_EXPRESSION(sx, sx_str);
    INIT_EXPRESSION(sy, sy_str);
    INIT_EXPRESSION(gx, gx_str);
    INIT_EXPRESSION(gy, gy_str);
}

void segy_free_reader_par(segy_reader_params_t* p)
{
    if(p->sx_str)             free(p->sx_str);
    if(p->sx) destroy_trchead_eval(p->sx);
    if(p->sy_str)             free(p->sy_str);
    if(p->sy) destroy_trchead_eval(p->sy);
    if(p->gx_str)             free(p->gx_str);
    if(p->gx) destroy_trchead_eval(p->gx);
    if(p->gy_str)             free(p->gy_str);
    if(p->gy) destroy_trchead_eval(p->gy);

    if(p->sepweights) sep_close(p->sepweights);
    if(p->weights) free(p->weights);

    if(p->mmap_flags) free(p->mmap_flags);

    memset(p, 0, sizeof(*p));
}

static char* segy_get_next_line_in_file_list(FILE* f)
{
    char* line;
    while((line = se_read_text_stream_line(f, 512)) != NULL) {
        char* ptr = strchr(line, '#');
        size_t len;

        if(ptr) *ptr = 0;
        str_trim(line);
        len = strlen(line);
        if( len > 0 ) return line;
        free(line);
    }
    return NULL;
}

static char** segy_get_file_list(const char* file_name)
{
    char** ret = NULL;

    if(file_name && se_file_exists(file_name)) {
        int is_text_list = 0;
        FILE *f = se_fopen(file_name, "r");
        if(f) {
            char* line = segy_get_next_line_in_file_list(f);
            if(line) {
                if( se_file_exists(line) ) {
                    array list = create_array();
                    is_text_list = 1;
                    array_add(list, line);
                    while( (line = segy_get_next_line_in_file_list(f)) != NULL ) {
                        array_add(list, line);
                    }
                    ret = array_to_null_ended_str_list(list);
                    destroy_array(list, 0);
                } else {
                    free(line);
                }
            }
            fclose(f);
        }

        if(!is_text_list) {
            int count = 1;
            array list = create_array();
            char* fn = asprintf("%s%d", file_name, count);

            array_add(list, (void *)strdup(file_name));
            while(se_file_exists(fn)) {
                array_add(list, fn);
                fn = asprintf("%s%d", file_name, ++count);
            }
            free(fn);

            ret = array_to_null_ended_str_list(list);

            destroy_array(list, 0);
        }
    }

    if(!ret)
    {
        ret = (char**)malloc(sizeof(*ret));
        memset(ret, 0, sizeof(*ret));
    }

    return ret;
}

static segy_file* segy_open_data_file(const extract_traces* et, const char* file_name)
{
    binhead bhead;
    trchead thead;
    segy_file *sf;
    int64_t fsize;

    /* Check if the current data file exists */
    if( ! se_file_exists(file_name) ) {
        ERROR(( "Cannot read from input file %s", file_name ));
        return NULL;
    }

    sf = (segy_file*)malloc(sizeof(*sf));
    memset(sf, 0, sizeof(*sf));
    sf->file_name = se_strdup(file_name);

    INFOV((11, "Opening data file %s", sf->file_name));
    sf->io_data = se_fsio_init(sf->file_name, "r", NULL);

    se_fsio_raw_read(sf->io_data, sf->reel, EBC_BYTES+BHEAD_BYTES);
    memcpy(&bhead, sf->reel+EBC_BYTES, BHEAD_BYTES);
    segy_binheader_order_bytes(&bhead);

    if( et->rp.forced_format == 0 ) sf->format = bhead.format;
    else                            sf->format = et->rp.forced_format;

    switch(sf->format) {
    case 1:
        sf->bps = 4;
        sf->fix_data = segy_fix_data_ibm_float;
        sf->ieee2segy = segy_ieee2segy_ibm_float;
        break;
    case 2:
        sf->bps = 4;
        sf->fix_data = segy_fix_data_fixed_point_4B;
        sf->ieee2segy = segy_ieee2segy_fixed_point_4B;
        break;
    case 3:
        sf->bps = 2;
        sf->fix_data = segy_fix_data_fixed_point_2B;
        sf->ieee2segy = segy_ieee2segy_fixed_point_2B;
        break;
    case 4:
        sf->bps = 4;
        sf->fix_data = segy_fix_data_fixed_point_gain;
        sf->ieee2segy = segy_ieee2segy_fixed_point_gain;
        break;
    case 5:
        sf->bps = 4;
        sf->fix_data = segy_fix_data_ieee_float;
        sf->ieee2segy = segy_ieee2segy_ieee_float;
        break;
    case 6:
        sf->bps = 1;
        sf->fix_data = segy_fix_data_fixed_point_1B;
        sf->ieee2segy = segy_ieee2segy_fixed_point_1B;
        break;
    default:
        WARN(( "Unknown data format code %d in file %s. Using 1.", sf->format, sf->file_name));
        sf->format = 1;
        sf->bps = 4;
        sf->fix_data = segy_fix_data_ibm_float;
        sf->ieee2segy = segy_ieee2segy_ibm_float;
    }

    se_fsio_raw_read(sf->io_data, &thead, THEAD_BYTES);
    segy_trcheader_order_bytes(&thead);

    if( thead.ns != bhead.hns ) {
        WARN(( "File %s: number of samples read from binary header (%d) "
               "is not the same as the number of samples read from first "
               "trace header (%d).",
               sf->file_name, bhead.hns, thead.ns ));
    }
    if(et->rp.forced_ns) {
        sf->ns = et->rp.forced_ns;
    } else {
        sf->ns = thead.ns;
    }

    if(et->rp.forced_idt) {
        sf->idt = et->rp.forced_idt;
    } else {
        sf->idt = thead.dt;
    }

    fsize = se_file_length(sf->file_name);
    sf->ntraces = (int)( (fsize - EBC_BYTES - BHEAD_BYTES) /
                         ( THEAD_BYTES + sf->ns * sf->bps ) );

    if( (int64_t)sf->ntraces * ( THEAD_BYTES + sf->ns * sf->bps ) !=
        (fsize - EBC_BYTES - BHEAD_BYTES) ) {
        WARN(("File %s doesn't have an exact number of traces", sf->file_name));
    }

    sf->trace = alloc1type(byte, THEAD_BYTES + sf->ns * sf->bps );
    sf->trace_data = alloc1float( sf->ns );

    INFOV((100, "Found file %s with %d traces, ns=%d, dt=%d",
           sf->file_name, sf->ntraces, sf->ns, sf->idt));

    se_fsio_close(sf->io_data); sf->io_data = NULL;

    return sf;
}

static extract_traces* init_segy_read(const char* file_name,
                                      const segy_reader_params_t* rp,
                                      const char* parset)
{
    extract_traces* et;
    char** file_list;
    int i;

    et = (extract_traces*)malloc(sizeof(extract_traces));
    memset(et, 0, sizeof(*et));
    et->last_sf_open = NULL;

    if(rp) {
        et->rp = *rp;
        segy_init_reader_par_expression_strings(&(et->rp),
                                                   rp->sx_str,
                                                   rp->sy_str,
                                                   rp->gx_str,
                                                   rp->gy_str);
    } else {
        segy_init_reader_par_from_parset(&(et->rp), parset);
    }


    if(et->rp.forced_format) {
        INFOV((10, "Forcing format code to %d for %s", et->rp.forced_format, file_name));
    }

    if(se_get_verb_level() >= 10) {
        switch( et->rp.use_scalco ) {
        case 0:
            INFO(("Ignoring scale factor. Using directly coordinates from %s", file_name));
            break;
        case 1:
            INFO(("Scale factor from scalco field, used as a multiplier for coordinates from %s", file_name));
            break;
        case 2:
            INFO(("Scale factor from scalco field, used as a power of 10 for coordinates from %s", file_name));
            break;
        default:
            INFO(("Using scale factor %f for coordinates from %s", et->rp.scalco_factor, file_name));
            break;
        }
    }

    if(et->rp.forced_ns)  {
        INFOV((10, "Forcing number of samples to %d for %s",
               et->rp.forced_ns, file_name));
        if(et->rp.forced_ns < 0) {
            ERROR(("Number of samples was set to a negative number: %d",
                   et->rp.forced_ns));
        }
    }

    if(et->rp.ns)  {
        INFOV((10, "Using %d samples from traces for %s",
               et->rp.ns, file_name));
        if(et->rp.ns < 0) {
            ERROR(("Number of samples to use was set to a negative number: %d",
                   et->rp.ns));
        }
    }

    if(et->rp.forced_idt) {
        INFOV((10, "Forcing dt to %d for %s", et->rp.forced_idt, file_name));
    }

    if(et->rp.sx_str && et->rp.sx)
        INFOV((10, "Replacing sx with [%s] when reading from %s", et->rp.sx_str, file_name));
    if(et->rp.sy_str && et->rp.sy)
        INFOV((10, "Replacing sy with [%s] when reading from %s", et->rp.sy_str, file_name));
    if(et->rp.gx_str && et->rp.gx)
        INFOV((10, "Replacing gx with [%s] when reading from %s", et->rp.gx_str, file_name));
    if(et->rp.gy_str && et->rp.gy)
        INFOV((10, "Replacing gy with [%s] when reading from %s", et->rp.gy_str, file_name));


    file_list = segy_get_file_list(file_name);


    for(et->nfiles = 0; file_list[et->nfiles]; et->nfiles += 1)
        ; // 空操作，只计算文件数量

    et->data = (segy_file**)malloc(et->nfiles*sizeof(*(et->data)));
    for(i = 0; i < et->nfiles; ++i) {
        et->data[i] = segy_open_data_file(et, file_list[i]);
        if(i > 0 ) {
            if(et->data[i]->ns != et->data[i-1]->ns) {
                ERROR(("Number of samples on file %s (%d) is different than for "
                       "the previous one %s (%d)",
                       et->data[i]->file_name  , et->data[i]->ns,
                       et->data[i-1]->file_name, et->data[i-1]->ns));
            }
            if( et->data[i]->idt != et->data[i-1]->idt ) {
                ERROR(( "Sampling on file %s (%d) is different than for "
                        "the previous one %s (%d)",
                        et->data[i]->file_name  , et->data[i]->idt,
                        et->data[i-1]->file_name, et->data[i-1]->idt));
            }
        }
    }

    free_null_ended_str_list(file_list);

    if( et->nfiles == 0 ) {
        ERROR(("Cannot open data file(s): %s", file_name));
    }

    if(et->rp.ns == 0) {
        et->rp.ns = et->data[0]->ns;
    }

    if(et->rp.weights) {
        int ntraces = segy_get_ntraces((segy_read_t)et);
        et->rp.sepweights = sep_open(et->rp.weights, SEP_READ, 0);
        if(et->rp.sepweights->headers->n[1] != ntraces) {
            ERROR(("Unexepectd number of traces (n1=%d) in the weights file %s; it should be %d",
                   et->rp.sepweights->headers->n[0], et->rp.weights, ntraces));
        }
    }

    return et;
}



void segy_destroy_segy_reader(segy_read_t segy)
{
    extract_traces* et = (extract_traces*)segy;
    int i;
    for(i = 0; i < et->nfiles; i++) {
        free(et->data[i]->file_name);
        free(et->data[i]->trace);
        free(et->data[i]->trace_data);
        if(et->data[i]->io_data) se_fsio_close(et->data[i]->io_data);
        free(et->data[i]);
    }
    if(et->nfiles) free(et->data);

    segy_free_reader_par(&(et->rp));

    free(et);
}




void segy_compute_trace_coordinates( const segy_read_t segy,
                                        const trchead* head,
                                        float* lsx, float* lsy,
                                        float* lgx, float* lgy )
{
    const extract_traces* et = (const extract_traces*)segy;
    double scale;
    int scalco = head->scalco;

    switch( et->rp.use_scalco ) {
    case 0:
        scale = 1.0;
        break;
    case 1:
        if(scalco > 1)       scale = scalco;
        else if(scalco < -1) scale = -1.0/(double)scalco;
        else                 scale = 1.0;
        break;
    case 2:
        scale = (double)pow( (double)10.0, scalco );
        break;
    default:
        scale = et->rp.scalco_factor;
        break;
    }

    if(et->rp.sx) *lsx = (float)trchead_eval(et->rp.sx, head);
    else          *lsx = (float)(scale*head->sx);
    if(et->rp.sy) *lsy = (float)trchead_eval(et->rp.sy, head);
    else          *lsy = (float)(scale*head->sy);
    if(et->rp.gx) *lgx = (float)trchead_eval(et->rp.gx, head);
    else          *lgx = (float)(scale*head->gx);
    if(et->rp.gy) *lgy = (float)trchead_eval(et->rp.gy, head);
    else          *lgy = (float)(scale*head->gy);
}

void segy_coords_to_integer(const segy_read_t segy,
                               const trchead* head,
                               double sx, double sy, double gx, double gy,
                               int* isx, int* isy, int* igx, int* igy)
{
    const extract_traces* et = (const extract_traces*)segy;
    double scale;
    int scalco = head->scalco;

    switch( et->rp.use_scalco ) {
    case 0:
        scale = 1.0;
        break;
    case 1:
        if(scalco > 1)       scale = 1.0/scalco;
        else if(scalco < -1) scale = -1.0 * (double)scalco;
        else                 scale = 1.0;
        break;
    case 2:
        scale = (double)pow( (double)10.0, -scalco );
        break;
    default:
        scale = 1.0/et->rp.scalco_factor;
        break;
    }
    *isx = roundi(scale*sx);
    *isy = roundi(scale*sy);
    *igx = roundi(scale*gx);
    *igy = roundi(scale*gy);
}

void segy_ieee2segy(const segy_read_t segy,
                       trchead* head,
                       byte* bdata, const float* data,
                       int nsamples)
{
    const extract_traces* et = (const extract_traces*)segy;
    const segy_file *sf = et->data[0];
    sf->ieee2segy(head, bdata, data, nsamples);
}

int segy_detect_null_trace(const float* tr, int64_t n){
	int64_t j;
	int ret=1;
	for(j=0; j < n; j++) if (!FISZERO(tr[j])) {
		ret=0;
		break;
	}
    return ret;
}
