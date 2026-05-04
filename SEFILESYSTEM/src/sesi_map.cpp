#include "../include/sesi_map.h"



#define SESI_COMPUTE_HASH_XOR

typedef struct _sesi_cell_cache_s _sesi_cell_cache_t;             //s structure t typedef  _sesi_cell_cache_t

typedef struct _sesi_map_s {
    int read_only;
    int naxis;
    axa_t* axis;
    int include_all;
    int grow_by;
    sep_t* sep;  //sepsepsepspesepssseepp

    int use_cache;
    size_t hashsize;
    int64_t max_cache_size;
    int64_t cache_size;
    _sesi_cell_cache_t** cache;
    _sesi_cell_cache_t* oldest;
    _sesi_cell_cache_t* newest;


    int64_t cache_lookups;
    int64_t cache_chain_walks;
    int64_t cache_n_elements;
    int64_t cache_n_inserts;
    int64_t cache_n_deletes;
    int64_t cache_n_read_info_hits;
    int64_t cache_n_read_info_misses;
    int64_t cache_n_write_info_hits;
    int64_t cache_n_write_info_misses;
    int64_t cache_n_read_list_hits;
    int64_t cache_n_read_list_misses;
    int64_t cache_n_write_list_hits;
    int64_t cache_n_write_list_misses;

} _sesi_map_t;

typedef struct _sesi_cell_s {
    /**
     * Stores the fold: how many indexes are in this cell.  If fold is
     * 0, 1 or 2, then the indexes are written to disk in the map info
     * area, packed to the place where the offset is otherwise: the
     * next two int32 entries.
     *
     * If the fold is bigger than 2, then the next two int32 entries
     * from the info map represents a 64-bit integer interpreted as
     * the offset where the actual list is.
     */
    int fold;
    
    /**
     * The offset in the file where this cell stores the list of
     * indexes. This is used only if the fold is bigger than 2,
     * meaning that the list is stored somewhere after the end of the
     * map info area.
     * 
     * Since the map info area is allways bigger than zero, a value of
     * zero for offset signifies that this cell doesn't have yet a
     * block allocated on disk.
     */
    int64_t off;

    /**
     * Used only in memory to as the list of indexes. When read from
     * this, here is stored the result of decoding the last two int32
     * from the info.
     */
    int idx[2];
}_sesi_cell_t;

struct _sesi_cell_cache_s {
    int64_t lidx;
    _sesi_cell_t cell;

    int list_size;
    int* list;

    _sesi_cell_cache_t* next;          //链表
    _sesi_cell_cache_t* newer;
};

static void sesi_map_discard_cache(_sesi_map_t* zsm);
static void sesi_map_flush_cache(_sesi_map_t* zsm);
static void sesi_map_flush_cell_cache(_sesi_map_t* zsm, _sesi_cell_cache_t* p);


sesi_map_t sesi_map_create(const char* file, const axa_t* axis, int naxis,
                             int include_all, int grow_by)
{
    int i;
     _sesi_map_t* zsm = (_sesi_map_t*)malloc(sizeof(*zsm));
    memset(zsm, 0, sizeof(*zsm));
    off_t map_size;

    zsm->sep = sep_open( file, SEP_WRITE, 0 );
    zsm->naxis = naxis;
    zsm->axis = (axa_t *)malloc(zsm->naxis * sizeof(axa_t));
    memset(zsm->axis, 0, zsm->naxis * sizeof(axa_t));
    for(i = 0; i < zsm->naxis; ++i) {
        zsm->axis[i] = axis[i];
    }

    zsm->sep->headers->ndim = naxis+1;
    zsm->sep->headers->o[0] = 0;
    zsm->sep->headers->d[0] = 1;
    zsm->sep->headers->n[0] = 3;
    map_size = 3*sizeof(float);
    for(i = 0; i < zsm->naxis; ++i) {
        zsm->sep->headers->o[i+1] = zsm->axis[i].o;
        zsm->sep->headers->d[i+1] = zsm->axis[i].d;
        zsm->sep->headers->n[i+1] = zsm->axis[i].n;
        map_size *= zsm->axis[i].n;
    }

    zsm->include_all = include_all;
    if(grow_by < 0) {
        grow_by = 31;
    }
    zsm->grow_by = grow_by;

    sep_set_header_int(zsm->sep, "sesi_map_include_all", zsm->include_all);
    sep_set_header_int(zsm->sep, "sesi_map_grow_by", zsm->grow_by);

    se_fsio_zero_file(zsm->sep->data->io, map_size);

    sep_write_headers( zsm->sep );


    zsm->read_only = 0;

    return (sesi_map_t)zsm;
}

sesi_map_t sesi_map_open(const char* file, int read_only)
{
    _sesi_map_t* zsm = (_sesi_map_t*)malloc(sizeof(*zsm));
    memset(zsm, 0, sizeof(*zsm));
    int i, sep_flag;
    off_t min_map_size;

    zsm->read_only = read_only;

    sep_flag = SEP_READ;
    if(!read_only) {
        sep_flag |= SEP_WRITE;
    }

    zsm->sep = sep_open( file, sep_flag, 0 );

    if(zsm->sep->headers->ndim < 2) {
        ERROR(("%s doesn't seem to be a SEP-map file (number of dimensions should be > 1)",
               file));
    }
    if(zsm->sep->headers->n[0] != 3) {
        ERROR(("%s doesn't seem to be a SEP-map file (n1 must be 3)",
               file));
    }

    zsm->include_all = sep_get_hdr_int(zsm->sep, "sesi_map_include_all", -1);
    if(zsm->include_all < 0) {
        ERROR(("%s may be a corrupted SEP-map file (wrong sesi_map_include_all)",
               file));
    }

    zsm->grow_by = sep_get_hdr_int(zsm->sep, "sesi_map_grow_by", -1);
    if(zsm->grow_by <= 0) {
        ERROR(("%s may be a corrupted SEP-map file (wrong sesi_map_grow_by)",
               file));
    }

    zsm->naxis = zsm->sep->headers->ndim-1;
    zsm->axis = (axa_t *)malloc(zsm->naxis * sizeof(axa_t));
    memset(zsm->axis, 0, zsm->naxis * sizeof(axa_t));
    min_map_size = 3*sizeof(float);
    for(i = 0; i < zsm->naxis; ++i) {
        zsm->axis[i].o = zsm->sep->headers->o[i+1];
        zsm->axis[i].d = zsm->sep->headers->d[i+1];
        zsm->axis[i].n = zsm->sep->headers->n[i+1];
        min_map_size *= zsm->axis[i].n;
    }

    if(min_map_size > se_fsio_length(zsm->sep->data->io)) {
        ERROR(("%s may be a corrupted SEP-map file (binary too small)",
               file));
    }

    return (sesi_map_t)zsm;
}

void sesi_map_print_cache_stat(sesi_map_t _zsm, int v)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    size_t i;
    int64_t nelem, nempty, max, min;
    se_quantiler_t zqt;

    if(!zsm->use_cache || se_get_verb_level() < v) return;

    INFO(("SI-MAP Cache Statistics:"));
    INFO(("\tcache_lookups=%Ld", (long long int)zsm->cache_lookups));
    INFO(("\tcache_chain_walks=%Ld", (long long int)zsm->cache_chain_walks));
    INFO(("\tcache_n_elements=%Ld", (long long int)zsm->cache_n_elements));
    INFO(("\tcache_n_inserts=%Ld", (long long int)zsm->cache_n_inserts));
    INFO(("\tcache_n_deletes=%Ld", (long long int)zsm->cache_n_deletes));
    INFO(("\tcache_n_read_info_hits=%Ld", (long long int)zsm->cache_n_read_info_hits));
    INFO(("\tcache_n_read_info_misses=%Ld", (long long int)zsm->cache_n_read_info_misses));
    INFO(("\tcache_n_write_info_hits=%Ld", (long long int)zsm->cache_n_write_info_hits));
    INFO(("\tcache_n_write_info_misses=%Ld", (long long int)zsm->cache_n_write_info_misses));
    INFO(("\tcache_n_read_list_hits=%Ld", (long long int)zsm->cache_n_read_list_hits));
    INFO(("\tcache_n_read_list_misses=%Ld", (long long int)zsm->cache_n_read_list_misses));
    INFO(("\tcache_n_write_list_hits=%Ld", (long long int)zsm->cache_n_write_list_hits));
    INFO(("\tcache_n_write_list_misses=%Ld", (long long int)zsm->cache_n_write_list_misses));

    nelem = 0;
    max = 0;
    min = 0;
    nempty = 0;

    se_quantiler_init(&zqt, 0.5, 0, 0);

    for(i = 0; i < zsm->hashsize; ++i) {
        _sesi_cell_cache_t* p = zsm->cache[i];
        int n = 0;
        while(p) {
            p = p->next;
            n += 1;
            nelem += 1;
        }
        if(n) {
            if(min) {
                if(max < n) max = n;
                if(min > n) min = n;
            } else {
                min = max = n;
            }
            se_quantiler_update(&zqt, (real)n);
        } else {
            nempty += 1;
        }
    }

    INFO(("\tcomputed stat: nelem=%Ld, min=%Ld, max=%Ld, hash=%Ld, empty=%Ld (%.2f%%) - median=%g",
          (long long int)nelem, (long long int)min, (long long int)max,
          (long long int)zsm->hashsize, (long long int)nempty, nempty*100.0/zsm->hashsize,
          se_quantiler_estimate(&zqt)));
}

int sesi_map_get_naxis(sesi_map_t _zsm)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    return zsm->naxis;
}

int64_t sesi_map_get_nelements(sesi_map_t _zsm)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    return sep_get_total_size(zsm->sep)/3;
}

axa_t sesi_map_get_axis(sesi_map_t _zsm, int i)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    ASSERT(i < zsm->naxis);
    return zsm->axis[i];
}

int sesi_map_have_hdr_int(sesi_map_t _zsm, const char* name)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    return sep_have_hdr_int(zsm->sep, name);
}
int sesi_map_get_hdr_int(sesi_map_t _zsm, const char* name, int def)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    return sep_get_hdr_int(zsm->sep, name, def);
}
int sesi_map_have_hdr_float(sesi_map_t _zsm, const char* name)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    return sep_have_hdr_float(zsm->sep, name);
}
double sesi_map_get_hdr_float(sesi_map_t _zsm, const char* name, double def)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    return sep_get_hdr_float(zsm->sep, name, def);
}
int sesi_map_have_hdr(sesi_map_t _zsm, const char* name)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    return sep_have_hdr(zsm->sep, name);
}
char* sesi_map_get_hdr(sesi_map_t _zsm, const char* name, const char* def)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    return sep_get_hdr(zsm->sep, name, def);
}

void sesi_map_set_header(sesi_map_t _zsm, const char* name, const char* value)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    sep_set_header(zsm->sep, name, value);
}
void sesi_map_set_header_int(sesi_map_t _zsm, const char* name, int value)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    sep_set_header_int(zsm->sep, name, value);
}
void sesi_map_set_header_float(sesi_map_t _zsm, const char* name, double value)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    sep_set_header_float(zsm->sep, name, value);
}

void sesi_map_set_header_survey_description(sesi_map_t _zsm,
                                             const survey_description_t* sd,
                                             const char* prefix,
                                             int use_only_set_fields)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    sep_set_header_survey_description(zsm->sep, sd, prefix, use_only_set_fields);
}

void sesi_map_get_header_survey_description(sesi_map_t _zsm,
                                             survey_description_t* sd,
                                             const char* prefix,
                                             const survey_description_t* project_sd)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    sep_get_header_survey_description(zsm->sep, sd, prefix, project_sd, 0);
}


void sesi_map_close(sesi_map_t _zsm)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;

    if(! zsm->read_only)
        sesi_map_flush_cache(zsm);

    sesi_map_discard_cache(zsm);
    if(zsm->cache) free(zsm->cache);

    free(zsm->axis);
    sep_close(zsm->sep);

    free(zsm);
}

inline int64_t sesi_map_get_cell_index(const _sesi_map_t* zsm, 
                                            const double* coords)
{
    int64_t lidx, nfac;
    int i;

    for(lidx = 0, nfac = 1, i = 0; i < zsm->naxis; ++i) {
        int64_t aidx;
        if(zsm->include_all) {
            aidx = axa_nidx_trim(zsm->axis+i, (real)coords[i]);
        } else {
            aidx = axa_nidx(zsm->axis+i, (real)coords[i]);
            if(aidx < 0 || aidx >= zsm->axis[i].n) return -1;
        }

        lidx += nfac*aidx;
        nfac *= zsm->axis[i].n;
    }

    return lidx;
}

inline int64_t sesi_map_get_number_of_cells(const _sesi_map_t* zsm)
{
    int i;
    int64_t ncell;

    for(ncell = 1, i = 0; i < zsm->naxis; ++i) {
        ncell *= zsm->axis[i].n;
    }
    return ncell;
}

int64_t sesi_map_get_number_of_non_empty_cells(sesi_map_t _zsm)
{
    int64_t i, n;
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    
    n=0;
    for (i=0; i<sesi_map_get_number_of_cells(zsm); i++) {
        if (sesi_map_get_fold_from_index(_zsm, i) > 0) 
            n++;
    }
    
    return n;
}


inline int64_t sesi_map_cell_cache_size(_sesi_cell_cache_t* p)
{
    return sizeof(*p) + p->list_size * sizeof(int);
}

inline size_t sesi_map_get_hash(int64_t lidx, size_t hashsize)
{
#ifdef SESI_COMPUTE_HASH_XOR
    return (size_t)( (
                      (lidx>>56) ^ (lidx>>48) ^ (lidx>>40) ^ (lidx>>32) ^
                      (lidx>>24) ^ (lidx>>16) ^ (lidx>> 8) ^ (lidx>> 0)
                      ) % hashsize );
#else
    return (size_t)(lidx % hashsize);
#endif
}

static _sesi_cell_cache_t* sesi_map_cache_lookup(_sesi_map_t* zsm, int64_t lidx)
{
    _sesi_cell_cache_t* p;
    zsm->cache_lookups += 1;
    for( p = zsm->cache[sesi_map_get_hash(lidx, zsm->hashsize)];
         p && p->lidx != lidx; p = p->next) zsm->cache_chain_walks += 1;
    return p;
}

static void sesi_map_cache_delete(_sesi_map_t* zsm, _sesi_cell_cache_t* p);

static void sesi_map_cache_destroy_cell_cache(_sesi_cell_cache_t *p)
{
    if(p->list) free(p->list);
    p->list = NULL;
    free(p);
}

static void sesi_map_cache_insert(_sesi_map_t* zsm, _sesi_cell_cache_t* p)
{
    size_t hash;

    ASSERT( p->next == 0 );

    hash = sesi_map_get_hash(p->lidx, zsm->hashsize);
    p->next = zsm->cache[hash];
    zsm->cache[hash] = p;

    p->newer = NULL;
    if(zsm->newest)
        zsm->newest->newer = p;
    zsm->newest = p;
    if(zsm->oldest == NULL) zsm->oldest = p;

    zsm->cache_size += sesi_map_cell_cache_size(p);
    while(zsm->cache_size > zsm->max_cache_size && zsm->oldest != p ) {
        _sesi_cell_cache_t *old = zsm->oldest;
        sesi_map_cache_delete(zsm, zsm->oldest);
        sesi_map_flush_cell_cache(zsm, old);
        sesi_map_cache_destroy_cell_cache(old);
    }

    zsm->cache_n_elements += 1;
    zsm->cache_n_inserts += 1;
}


static void sesi_map_cache_delete(_sesi_map_t* zsm, _sesi_cell_cache_t* p)
{
    _sesi_cell_cache_t **pp;

    pp = &zsm->cache[sesi_map_get_hash(p->lidx, zsm->hashsize)];
    for( ; (*pp) != p; pp = &(*pp)->next) { ASSERT(*pp); }
    *pp = p->next;
    p->next = 0;

    zsm->cache_size -= sesi_map_cell_cache_size(p);
    if(p == zsm->oldest) {
        zsm->oldest = p->newer;
    }

    zsm->cache_n_elements -= 1;
    zsm->cache_n_deletes += 1;
}


static void sesi_map_discard_cache(_sesi_map_t* zsm)
{
    if(zsm->use_cache) {
        size_t i;
        for(i = 0; i < zsm->hashsize; ++i) {
            _sesi_cell_cache_t* p = zsm->cache[i];
            while(p) {
                _sesi_cell_cache_t* old = p;
                p = p->next;
                zsm->cache_size -= sesi_map_cell_cache_size(old);
                sesi_map_cache_destroy_cell_cache(old);
            }
            zsm->cache[i] = NULL;
        }
        zsm->newest = zsm->oldest = NULL;
    }
}

/**
 * Returns the cell info structure for the lidx. If not found in
 * cache, the info is read and decoded from the map info area of the
 * file.
 *
 * If necessary, a new entry to the cache is added.
 */
static void sesi_map_read_cell_info(_sesi_map_t* zsm,
                                          se_fsio* io, int64_t lidx, _sesi_cell_t* cell,
                                          int use_cache, int create_new_cache)
{
    _sesi_cell_cache_t* pcache = NULL;
    union {
        int32_t i[3];
        float f[3];
    } buf;

    if(use_cache && zsm->use_cache) {
        pcache = sesi_map_cache_lookup(zsm, lidx);
    }

    if( pcache ) { /*if in cache, just copy the cell info from pcache */
        *cell = pcache->cell;
        zsm->cache_n_read_info_hits += 1;
    } else {
        se_fsio_seek(io, lidx*3*sizeof(float));
        se_fsio_read_float(io, buf.f, 3);

        cell->fold = (int)buf.f[0];
        if(cell->fold < 0) {
            ERROR(("Corrupted map file - negative fold found"));
        }

        if(cell->fold == 0) {
            cell->off = 0;
            cell->idx[0] = -1;
            cell->idx[1] = -1;
        } else if(cell->fold == 1) {
            cell->off = 0;
            cell->idx[0] = buf.i[1];
            cell->idx[1] = -1;
        } else if(cell->fold == 2) {
            cell->off = 0;
            cell->idx[0] = buf.i[1];
            cell->idx[1] = buf.i[2];
        } else {
            cell->off = (((int64_t)buf.i[1])<<32)|((uint32_t)buf.i[2]);
            cell->idx[0] = -1;
            cell->idx[1] = -1;
        }

        zsm->cache_n_read_info_misses += 1;
    }

    if(create_new_cache && use_cache && zsm->use_cache && !pcache) {
        pcache = (_sesi_cell_cache_t *)malloc(sizeof(*pcache));
        pcache->lidx = lidx;
        pcache->cell = *cell;
        pcache->list_size = 0;
        pcache->list = NULL;
        pcache->next = pcache->newer = NULL;
        sesi_map_cache_insert(zsm, pcache);
    }
}

/**
 * Writes the cell info to the map info area. If the fold is bigger
 * than 2, then the offset from the cell must be correctly specified
 * at this point.
 *
 * If the cell is already in cache, then the entry in cache is updated.
 */
static void sesi_map_write_cell_info(_sesi_map_t* zsm,
                                           se_fsio* io, int64_t lidx, const _sesi_cell_t* cell,
                                           int use_cache)
{
    _sesi_cell_cache_t* pcache = NULL;
    union {
        int32_t i[3];
        float f[3];
    } buf;

    if(use_cache && zsm->use_cache) {
        pcache = sesi_map_cache_lookup(zsm, lidx);
    }

    if( pcache ) { /* if in cache, update the entry */
        if(&(pcache->cell) != cell) pcache->cell = *cell;
        zsm->cache_n_write_info_hits += 1;
    } else {
        if(cell->fold == 0) {
            buf.f[0] = 0;
            buf.i[1] = -1;
            buf.i[2] = -1;
        } else if(cell->fold == 1) {
            buf.f[0] = 1;
            buf.i[1] = cell->idx[0];
            buf.i[2] = -1;
        } else if(cell->fold == 2) {
            buf.f[0] = 2;
            buf.i[1] = cell->idx[0];
            buf.i[2] = cell->idx[1];
        } else {
            buf.f[0] = (float)cell->fold;
            buf.i[1] = (int32_t)((cell->off)>>32);
            buf.i[2] = (uint32_t)(cell->off);
        }

        se_fsio_seek(io, lidx*3*sizeof(float));
        se_fsio_write_float(io, buf.f, 3);
        zsm->cache_n_write_info_misses += 1;
    }

    if(use_cache && zsm->use_cache && !pcache) {
        pcache = (_sesi_cell_cache_t *)malloc(sizeof(*pcache));
        pcache->lidx = lidx;
        pcache->cell = *cell;
        pcache->list_size = 0;
        pcache->list = NULL;
        pcache->next = pcache->newer = NULL;
        sesi_map_cache_insert(zsm, pcache);
    }
}


/**
 * Reads the list of indexes associated with the cell. If the cell is
 * already in cache, then the list is copied from there; otherwise the
 * list is read from this, at the location pointed by cell->off.
 *
 * The list is assumed to be allocated already and to be big enough to
 * hold fold+1 integers.
 */
static void sesi_map_read_cell_list(_sesi_map_t* zsm,
                                          se_fsio* io, 
                                          const _sesi_cell_t* cell, 
                                          int* list, int use_cache, int64_t lidx,
                                          int create_new_cache)
{
    _sesi_cell_cache_t* pcache = NULL;

    if(use_cache && zsm->use_cache) {
        pcache = sesi_map_cache_lookup(zsm, lidx);
    }

    if( pcache && pcache->list ) {
        ASSERT(pcache->list_size >= cell->fold + 1);
        memcpy(list, pcache->list, (cell->fold + 1)*sizeof(int));
        zsm->cache_n_read_list_hits += 1;
    } else {
        se_fsio_seek(io, cell->off);
        se_fsio_read_int32(io, list, cell->fold + 1);
        if(list[cell->fold] < cell->fold) {
            ERROR(("Corrupted map file: cell disk size is %d but fold is %d",
                   list[0], cell->fold));
        }
        zsm->cache_n_read_list_misses += 1;
    }

    if(use_cache && zsm->use_cache && (!pcache || !pcache->list)) {
        if(create_new_cache && !pcache) {
            pcache = (_sesi_cell_cache_t *)malloc(sizeof(*pcache));
            pcache->lidx = lidx;
            pcache->cell = *cell;
            pcache->list_size = 0;
            pcache->list = NULL;
            pcache->next = pcache->newer = NULL;
            sesi_map_cache_insert(zsm, pcache);
        }

        if(pcache) {
            if(pcache->list_size < pcache->cell.fold + 2) {
                zsm->cache_size -= pcache->list_size*sizeof(int);
                
                pcache->list_size = pcache->cell.fold + 2 + zsm->grow_by/2;
                int* new_list = realloc1int(pcache->list, pcache->list_size);
                if (new_list == NULL) {
                    ERROR(("Failed to reallocate cache list"));
                    return;
                }
                pcache->list = new_list;

                zsm->cache_size += pcache->list_size*sizeof(int);
            }
            
            memcpy(pcache->list, list, (cell->fold + 1)*sizeof(int));
        }
    }
}

/** 
 * Saves the cell info and list. If the cell is in cache, it just
 * updates the cache entry; otherwise the list is writen to disk.
 *
 * If write_extra si true, then when writing to the file, it also
 * writes zeros the the extra part of the block, between the end of
 * the list (fold+1) and the end of the block.
 */
static void sesi_map_write_cell_list(_sesi_map_t* zsm,
                                           se_fsio* io, 
                                           const _sesi_cell_t* cell, 
                                           int* list, int write_extra,
                                           int use_cache, int64_t lidx)
{
    _sesi_cell_cache_t* pcache = NULL;

    if(use_cache && zsm->use_cache) {
        pcache = sesi_map_cache_lookup(zsm, lidx);
    }

    /* if in cache, update the cache entry with the new info; adjust
       the memory for the list if needed.*/
    if( pcache ) {
        pcache->cell = *cell;
        if(pcache->list_size < pcache->cell.fold + 2) {
            zsm->cache_size -= pcache->list_size*sizeof(int);

            pcache->list_size = pcache->cell.fold + 2 + zsm->grow_by/2;
            int* new_list = realloc1int(pcache->list, pcache->list_size);
            if (new_list == NULL) {
                ERROR(("Failed to reallocate cache list in sesi_map"));
                return;
            }
            pcache->list = new_list;

            zsm->cache_size += pcache->list_size*sizeof(int);
        }

        memcpy(pcache->list, list, (cell->fold + 1)*sizeof(int));
        zsm->cache_n_write_list_hits += 1;
    } else {
        se_fsio_seek(io, cell->off);
        se_fsio_write_int32(io, list, cell->fold + 1);
        if(write_extra && list[cell->fold] > cell->fold) {
            int* extra = alloc1int_zero(list[cell->fold] - cell->fold);
            se_fsio_write_int32(io, extra, list[cell->fold] - cell->fold);
            free(extra);
        }
        zsm->cache_n_write_list_misses += 1;
    }
}

/**
 * Writes the information from the cache entry associated with a cell
 * to the disk. This function is called when flushing the entire cache
 * or when an old cache entry is selected to be discarded.
 */
static void sesi_map_flush_cell_cache(_sesi_map_t* zsm, _sesi_cell_cache_t* p)
{
    if(!zsm->use_cache || zsm->read_only)
        return; /*nothing to do*/

    if(p->cell.fold < 3) {
        ASSERT(p->cell.off == 0);
        sesi_map_write_cell_info(zsm, zsm->sep->data->io, 
                                  p->lidx, &p->cell, 0);
    } else {
        ASSERT(p->list); /* at this point, we need to have a list: fold bigger than 2 */

        p->cell.off = se_fsio_length(zsm->sep->data->io); //todo:review
        ASSERT(p->list_size > p->cell.fold);
        p->list[p->cell.fold] = p->cell.fold;

        sesi_map_write_cell_list(zsm, zsm->sep->data->io,
                                  &p->cell, 
                                  p->list, 0, 0, p->lidx);

        sesi_map_write_cell_info(zsm, zsm->sep->data->io, 
                                  p->lidx, &p->cell, 0);
    }
}

/**
 * Flushes the entire cache: writes to the disk all cell entries
 * currently held in cache.
 */
static void sesi_map_flush_cache(_sesi_map_t* zsm)
{
    if(zsm->use_cache) {
        size_t i;
        INFOV((100, "Flushing %Ld bytes of cache", (long long int)zsm->cache_size));
        for(i = 0; i < zsm->hashsize; ++i) {
            _sesi_cell_cache_t* p = zsm->cache[i];
            while(p) {
                sesi_map_flush_cell_cache(zsm, p);
                p = p->next;
            }
        }
    }
}

/**
 * Updates the list of ids for the lidx cell. If the append flag is
 * on, then the new list is added to the existing one. Othewise the
 * function replaces the content of the list.
 */
static int64_t sesi_map_update_many_internal (_sesi_map_t* zsm, int64_t lidx,
                                                    const int* idx, int nidx,
                                                    int append)
{
    _sesi_cell_t cell;
    int* list;
    int write_extra, block_size;
    int i;
    int start_fold;

    /* if no idexes are to be added, return now */
    if ( (append) && (nidx <= 0)) return lidx;

    /* read the cell info at lidx, and add it to the cache is not already there */
    sesi_map_read_cell_info(zsm, zsm->sep->data->io, lidx, &cell, zsm->use_cache, 1);

    start_fold = append ? cell.fold : 0;

    /* if the existing fold plus the new indexes gets to 0, 1 or 2*/
    if(start_fold + nidx < 3) {
        cell.fold = start_fold;
        for(i = 0; i < nidx; ++i) {
            cell.idx[cell.fold++] = idx[i];
        }

        /* just write the info: if the entry is in cache, then simply update info part of the entry */
        sesi_map_write_cell_info(zsm, zsm->sep->data->io, lidx, &cell, zsm->use_cache);
        return lidx;
    }

    write_extra = 0;
    list = alloc1int(start_fold + nidx + 1);


    /* first we copy the existing indexes into the list; the list may
       be either to the info block or in a list */
    if(start_fold < 3) { /* existing list is only in the cell info */
        for(i = 0; i < start_fold; ++i) {
            list[i] = cell.idx[i];
        }
        block_size = 0;
    } else {
        sesi_map_read_cell_list(zsm, zsm->sep->data->io, &cell, list, zsm->use_cache, lidx, 1);
        block_size = list[cell.fold];
    }

    cell.fold = start_fold;
    for(i = 0; i < nidx; ++i) {
        list[cell.fold++] = idx[i];
    }

    if(block_size < cell.fold) {
        block_size += (1 + (cell.fold - block_size + 1)/zsm->grow_by)*zsm->grow_by;
        cell.off = se_fsio_length(zsm->sep->data->io);
        write_extra = 1;
    }

    list[cell.fold] = block_size;

    sesi_map_write_cell_list(zsm, zsm->sep->data->io, &cell, list, write_extra, zsm->use_cache, lidx);
    sesi_map_write_cell_info(zsm, zsm->sep->data->io, lidx, &cell, zsm->use_cache);

    free(list);
    return lidx;
}

static int64_t sesi_map_add_internal (_sesi_map_t* zsm, int idx, int64_t lidx)
{
    return sesi_map_update_many_internal(zsm, lidx, &idx, 1, 1);
}


int64_t sesi_map_add_from_coords(sesi_map_t _zsm, int idx, const double* coords)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    int64_t lidx = sesi_map_get_cell_index(zsm, coords);
    if(lidx < 0) return -1;
    return sesi_map_add_internal(zsm, idx, lidx);
}

int64_t sesi_map_add_from_index(sesi_map_t _zsm, int idx, int64_t lidx)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    if( lidx < 0 || lidx >= sesi_map_get_number_of_cells(zsm) ) 
        return -1;
    return sesi_map_add_internal(zsm, idx, lidx);
}

int64_t sesi_map_set_indexes_from_coords(sesi_map_t _zsm, const double* coords,
                                          const int* idx, int nidx)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    int64_t lidx = sesi_map_get_cell_index(zsm, coords);
    if(lidx < 0) return -1;
    return sesi_map_update_many_internal(zsm, lidx, idx, nidx, 0);
}

int64_t sesi_map_set_indexes_from_index(sesi_map_t _zsm, int64_t lidx,
                                         const int* idx, int nidx)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    if( lidx < 0 || lidx >= sesi_map_get_number_of_cells(zsm) )
        return -1;
    return sesi_map_update_many_internal(zsm, lidx, idx, nidx, 0);
}

int64_t sesi_map_add_indexes_from_index(sesi_map_t _zsm,int64_t lidx,
                                        const int* idx, int nidx)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    if( lidx < 0 || lidx >= sesi_map_get_number_of_cells(zsm) )
        return -1;
    return sesi_map_update_many_internal(zsm, lidx, idx, nidx, 1);
}


void sesi_map_compact(sesi_map_t _zsm)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    int64_t lidx, ncells;
    char* tmp_filename;
    se_fsio* tmpio;
    int tmpfd;

    ncells = sesi_map_get_number_of_cells(zsm);

    tmp_filename = asprintf("%s.temp.XXXXXX", zsm->sep->data->filename);
    tmpfd = mkstemp(tmp_filename);
    if(tmpfd < 0) {
        ERROR(("Cannot open temp file %s", tmp_filename));
    }
    tmpio = se_fsio_attach( tmp_filename, tmpfd, zsm->sep->data->io->data_format);
    INFOV((100, "Initializing temp file"));
    se_fsio_zero_file(tmpio, ncells*3*sizeof(float));

    INFOV((100, "Start reading bins"));
    for(lidx = 0; lidx < ncells; ++lidx) {
        _sesi_cell_t cell;

        sesi_map_read_cell_info(zsm, zsm->sep->data->io, lidx, &cell, zsm->use_cache, 0);

        if(cell.fold > 2) {
            int* list = alloc1int(cell.fold + 1);
            sesi_map_read_cell_list(zsm, zsm->sep->data->io, &cell, list, zsm->use_cache, lidx, 0);
            list[cell.fold] = cell.fold;
            cell.off = se_fsio_length(tmpio);
            sesi_map_write_cell_list(zsm, tmpio, &cell, list, 0, 0, lidx);
            free(list);
        }

        sesi_map_write_cell_info(zsm, tmpio, lidx, &cell, 0);
    }
    INFOV((100, "Transfer temp file to final"));
    se_fsio_copy_file( tmpio, zsm->sep->data->io );
    se_fsio_remove(tmpio);
    free(tmp_filename);

    sesi_map_discard_cache(zsm);
}

int sesi_map_get_fold_from_coords(sesi_map_t _zsm, 
                                   const double* coords)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    _sesi_cell_t cell;
    int64_t lidx = sesi_map_get_cell_index(zsm, coords);
    if(lidx < 0) {
        return 0;
    }

    sesi_map_read_cell_info(zsm, zsm->sep->data->io, lidx, &cell, zsm->use_cache, 1);
    return cell.fold;
}

int sesi_map_get_fold_from_index(sesi_map_t _zsm, 
                                  int64_t lidx)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    _sesi_cell_t cell;
    if( lidx < 0 || lidx >= sesi_map_get_number_of_cells(zsm) ) {
        return 0;
    }

    sesi_map_read_cell_info(zsm, zsm->sep->data->io, lidx, &cell, zsm->use_cache, 1);
    return cell.fold;
}


static void sesi_map_get_indexes_from_index_internal(_sesi_map_t* zsm, 
                                                           int64_t lidx,
                                                           int** list, int* n, int* list_size)
{
    _sesi_cell_t cell;
    sesi_map_read_cell_info(zsm, zsm->sep->data->io, lidx, &cell, zsm->use_cache, 1);

    if(cell.fold == 0) {
        *n = 0;
        return;
    }

    if(*list_size < cell.fold + 2) {
        int* new_list = realloc1int(*list, cell.fold + 2);
        if (new_list == NULL) {
            ERROR(("Failed to reallocate index list"));
            return;
        }
        *list = new_list;
        *list_size = cell.fold + 2;
    }

    if(cell.fold == 1) {
        (*list)[0] = cell.idx[0];
        (*list)[1] = 0;
    } else if(cell.fold == 2) {
        (*list)[0] = cell.idx[0];
        (*list)[1] = cell.idx[1];
        (*list)[2] = 0;
    } else {
        sesi_map_read_cell_list(zsm, zsm->sep->data->io, &cell, *list, zsm->use_cache, lidx, 1);
    }

    *n = cell.fold;
}


void sesi_map_get_indexes_from_index(sesi_map_t _zsm, 
                                      int64_t lidx,
                                      int** list, int* n, int* list_size)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    if( lidx < 0 || lidx >= sesi_map_get_number_of_cells(zsm) ) {
        *n = 0;
        return;
    }

    sesi_map_get_indexes_from_index_internal(zsm, lidx, list, n, list_size);
}

void sesi_map_get_indexes_from_coords(sesi_map_t _zsm, 
                                       const double* coords,
                                       int** list, int* n, int* list_size)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;
    int64_t lidx = sesi_map_get_cell_index(zsm, coords);
    if(lidx < 0) {
        *n = 0;
        return;
    }

    sesi_map_get_indexes_from_index_internal(zsm, lidx, list, n, list_size);
}


void sesi_map_set_cache(sesi_map_t _zsm, int64_t max_cache_size, size_t hashsize)
{
    _sesi_map_t* zsm = (_sesi_map_t*)_zsm;

    if(zsm->use_cache) {
        if(! zsm->read_only)
            sesi_map_flush_cache(zsm);
        sesi_map_discard_cache(zsm);
        free(zsm->cache);
    }

    zsm->use_cache = 1;
    zsm->hashsize = hashsize;
    zsm->max_cache_size = max_cache_size;
    zsm->cache_size = 0;
    zsm->cache = (_sesi_cell_cache_t **)malloc(zsm->hashsize*sizeof(*zsm->cache));
    zsm->cache_size += zsm->hashsize*sizeof(*zsm->cache);
    zsm->newest = zsm->oldest = NULL;
}


int sesi_map_is_good_map_header(const char* file)
{
    return sep_is_good_sep_header(file);
}
