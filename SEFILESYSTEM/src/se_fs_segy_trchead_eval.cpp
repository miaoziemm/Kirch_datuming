#include "se_fs_segy.h"


/* 2*(240-4) 4b integers - for i* and ui* */
/* 2*(240-2) 2b integers - for s* and us* */
/* 2*(240-1) 1b integers - for b* and ub* */
/* 2*(240-4) 4b floating point - for fi* and f* */
/* 6 - for lsx,lsy,lgx,lgy,a,ra */
#define NEXTRAIDX (THEAD_NKEYS + 1898)
#define NVARS (NEXTRAIDX + 6)

typedef enum {
    OHT_INT,
    OHT_UINT,
    OHT_SHORT,
    OHT_USHORT,
    OHT_BYTE,
    OHT_UBYTE,
    OHT_FIBM,
    OHT_FIEEE
} offset_header_type;

typedef struct trc_eval_s {
    se_exp_t expr;
    double vars[NVARS];

    /* optimization of offset-only expressions */
    int off;
    offset_header_type oht;
} trc_eval;
static char* var_names[NVARS];
static int var_names_initialized = 0;
static char names_arena[ (240-4)*(1+3+1) + /* i*  */
                         (240-4)*(2+3+1) + /* ui* */
                         (240-2)*(1+3+1) + /* s*  */
                         (240-2)*(2+3+1) + /* us* */
                         (240-1)*(1+3+1) + /* b*  */
                         (240-1)*(2+3+1) + /* ub* */
                         (240-4)*(2+3+1) + /* fi* */
                         (240-4)*(1+3+1)   /* f* */
                         ];

#define HNAME_QUOTE(str) #str
#define IDX_NAME(fld) hdr_idx_##fld
#define SET_HDR_VAR_NAME(fld) var_names[hdr_idx_##fld]=strdup(#fld)

enum hdr_idx { SEGY_HEADER_FIELDS_LIST(IDX_NAME) };
/*static const char* hdr_names[] = { HEADER_FIELDS_LIST( HNAME_QUOTE ) };*/

static void init_offset_only_expression(trc_eval* te, const char* expr);
static double offset_only_evaluation(trc_eval* te, const byte* b);

trchead_eval_handle init_trchead_eval( const char* expr )
{
    int i;
    trc_eval* te;

    if( ! var_names_initialized ) {
        char* se_namespace;
        char** varname;

        for(i = 0; i < NVARS; ++i) var_names[i] = NULL;

        SEGY_HEADER_FIELDS_LIST(SET_HDR_VAR_NAME);

        se_namespace = names_arena;
        varname = &var_names[THEAD_NKEYS];
        for( i = 0; i < 240-4 ; i++, varname++, se_namespace += 1+3+1 ) {
            snprintf( *varname = se_namespace, 1+3+1,  "i%d", i+1 );
        }
        for( i = 0; i < 240-4 ; i++, varname++, se_namespace += 2+3+1 ) {
            snprintf( *varname = se_namespace, 2+3+1,  "ui%d", i+1 );
        }
        for( i = 0; i < 240-2 ; i++, varname++, se_namespace += 1+3+1 ) {
            snprintf( *varname = se_namespace, 1+3+1,  "s%d", i+1 );
        }
        for( i = 0; i < 240-2 ; i++, varname++, se_namespace += 2+3+1 ) {
            snprintf( *varname = se_namespace, 2+3+1,  "us%d", i+1 );
        }
        for( i = 0; i < 240-1 ; i++, varname++, se_namespace += 1+3+1 ) {
            snprintf( *varname = se_namespace, 1+3+1,  "b%d", i+1 );
        }
        for( i = 0; i < 240-1 ; i++, varname++, se_namespace += 2+3+1 ) {
            snprintf( *varname = se_namespace, 2+3+1,  "ub%d", i+1 );
        }
        for( i = 0; i < 240-4 ; i++, varname++, se_namespace += 2+3+1 ) {
            snprintf( *varname = se_namespace, 2+3+1,  "fi%d", i+1 );
        }
        for( i = 0; i < 240-4 ; i++, varname++, se_namespace += 1+3+1 ) {
            snprintf( *varname = se_namespace, 1+3+1,  "f%d", i+1 );
        }

        var_names_initialized = 1;

        ASSERT(varname - var_names == NEXTRAIDX);

        /*lsx,lsy,lgx,lgy,a,ra*/
        *(varname++) = strdup("lsx");
        *(varname++) = strdup("lsy");
        *(varname++) = strdup("lgx");
        *(varname++) = strdup("lgy");
        *(varname++) = strdup("a");
        *(varname++) = strdup("ra");

        ASSERT(varname - var_names == NVARS);

        for(i = 0; i < NVARS; ++i) {
            if(var_names[i] == NULL) ERROR(("Uninitialized SEGY eval at %d", i));
        }
    }


    if( expr == NULL ) {
        WARN(("Unexpected SEGY trace header expression: NULL"));
        return NULL;
    }

    te = (trc_eval*)malloc(sizeof(trc_eval));
    memset(te, 0, sizeof(*te));
    te->expr = se_exp_parse_simple( expr, (const char**)var_names, NVARS );
    if( se_exp_have_error(te->expr) ) {
        se_exp_print_formated_error(te->expr, 0, 0);
        free(te);
        return NULL;
    }

    init_offset_only_expression(te, expr);

    return ( trchead_eval_handle )te;
}

void destroy_trchead_eval( trchead_eval_handle eh )
{
    trc_eval* te = (trc_eval*)eh;

    if( te != NULL ) {
        if( te->expr != NULL )
            se_exp_destroy_expression(te->expr);

        free(te);
    }
}

void cleanup_trchead_eval_globals( void )
{
    if( var_names_initialized ) {
        int i;
        
        /* Clean up strdup allocated strings (only the ones after THEAD_NKEYS) */
        for( i = THEAD_NKEYS; i < NVARS; ++i ) {
            if( var_names[i] != NULL && i >= NEXTRAIDX ) {
                /* Only free the extra variables (lsx, lsy, lgx, lgy, a, ra) */
                free(var_names[i]);
                var_names[i] = NULL;
            }
        }
        
        /* Reset initialization flag */
        var_names_initialized = 0;
    }
}


#define B2I32(b1, b2, b3, b4) ( \
       ( ((uint32_t)(b1) << 24) & (uint32_t)0xff000000 ) | \
       ( ((uint32_t)(b2) << 16) & (uint32_t)0x00ff0000 ) | \
       ( ((uint32_t)(b3) <<  8) & (uint32_t)0x0000ff00 ) | \
       ( ((uint32_t)(b4)      ) & (uint32_t)0x000000ff ) )

#define B2I16(b1, b2) ( \
       ( ((uint16_t)(b1) <<  8) & (uint16_t)0xff00 ) | \
       ( ((uint16_t)(b2)      ) & (uint16_t)0x00ff ) )

#define GET_HDR_VAL(fld) te->vars[hdr_idx_##fld]=trc.fld

static void trchead_eval_fill_vars( trc_eval* te, const trchead* trch )
{
    trchead trc;
    int i;
    double* vars;

    /*memcpy( &trc, trch, sizeof(trchead) );*/
    trc = *trch;

    SEGY_HEADER_FIELDS_LIST(GET_HDR_VAL);

    segy_trcheader_order_bytes(&trc);

    vars = te->vars + THEAD_NKEYS;
    /* i* */
    for( i = 0; i < 240-4 ; i++, vars++ ) {
        byte* b = ( ((byte*)(&trc)) + i );
        int val = (int)B2I32(b[0], b[1], b[2], b[3]);
        *vars = (double)val;
    }
    /* ui* */
    for( i = 0; i < 240-4 ; i++, vars++ ) {
        byte* b = ( ((byte*)(&trc)) + i );
        unsigned int val = (unsigned int)B2I32(b[0], b[1], b[2], b[3]);
        *vars = (double)val;
    }
    /* s* */
    for( i = 0; i < 240-2 ; i++, vars++ ) {
        byte* b = ( ((byte*)(&trc)) + i );
        short val = (short)B2I16(b[0], b[1]);
        *vars = (double)val;
    }
    /* us* */
    for( i = 0; i < 240-2 ; i++, vars++ ) {
        byte* b = ( ((byte*)(&trc)) + i );
        unsigned short val = (unsigned short)B2I16(b[0], b[1]);
        *vars = (double)val;
    }
    /* b* */
    for( i = 0; i < 240-1 ; i++, vars++ ) {
        char val = *( (char*)( ((byte*)(&trc)) + i ) );
        *vars = (double)val;
    }
    /* ub* */
    for( i = 0; i < 240-1 ; i++, vars++ ) {
        unsigned char val = *( (unsigned char*)( ((byte*)(&trc)) + i ) );
        *vars = (double)val;
    }
    /* fi* */
    for( i = 0; i < 240-4 ; i++, vars++ ) {
        union {
            unsigned int i;
            float f;
        }i2f;
        unsigned int ival;
        byte* b = ( ((byte*)(&trc)) + i );
        ival = (unsigned int)B2I32(b[0], b[1], b[2], b[3]);
        ibm_to_ieee(&ival, &i2f.i, 1, 1);
        *vars = (double)i2f.f;
    }
    /* f* */
    for( i = 0; i < 240-4 ; i++, vars++ ) {
        byte* b = ( ((byte*)(&trc)) + i );
        union {
            unsigned int val;
            float f;
        }i2f;
        i2f.val = (unsigned int)B2I32(b[0], b[1], b[2], b[3]);
        *vars = (double)i2f.f;
    }

    for(i = NEXTRAIDX; i < NVARS; ++i) {
        te->vars[i] = 0;
    }
}

double trchead_eval( trchead_eval_handle eh, const trchead* trch )
{
    trc_eval* te = (trc_eval*)eh;

    if(te->off >= 0) {
        trchead trc = *trch;
        segy_trcheader_order_bytes(&trc);
        return offset_only_evaluation(te, (byte*)&trc);
    }


    trchead_eval_fill_vars(te, trch);
    return se_exp_eval_simle(te->expr, te->vars);
}

double trchead_evalx( trchead_eval_handle eh, const trchead* trch,
                          double lsx, double lsy, double lgx, double lgy)
{
    trc_eval* te = (trc_eval*)eh;
    int i;
    double a;

    if(te->off >= 0) {
        trchead trc = *trch;
        segy_trcheader_order_bytes(&trc);
        return offset_only_evaluation(te, (byte*)&trc);
    }

    trchead_eval_fill_vars(te, trch);

    a = segy_compute_azimuth(lsx, lsy, lgx, lgy);

    i = NEXTRAIDX-1;
    
    /*lsx,lsy,lgx,lgy,a,ra*/
    te->vars[++i] = lsx;
    te->vars[++i] = lsy;
    te->vars[++i] = lgx;
    te->vars[++i] = lgy;
    te->vars[++i] = a;
    te->vars[++i] = a*M_PI/180.0;

    return se_exp_eval_simle(te->expr, te->vars);
}

static void init_offset_only_expression(trc_eval* te, const char* expr)
{
    const char* c = expr;
    int u;
    int maxoff;

    te->off = -1;

    while(*c && isspace(*c)) c++;
    if(*c == 0) return;
    
    if(*c == 'u') {
        u = 1;
        c++;
    } else {
        u = 0;
    }

    if(*c == 0) return;

    switch(*c) {
    case 'b': 
        if(u) {
            te->oht = OHT_UBYTE;
        } else {
            te->oht = OHT_BYTE;
        }
        maxoff = THEAD_BYTES - 1;
        break;

    case 's':
        if(u) {
            te->oht = OHT_USHORT;
        } else {
            te->oht = OHT_SHORT;
        }
        maxoff = THEAD_BYTES - 2;
        break;

    case 'i': 
        if(u) {
            te->oht = OHT_UINT;
        } else {
            te->oht = OHT_INT;
        }
        maxoff = THEAD_BYTES - 4;
        break;

    case 'f':
        if(*(c+1) == 'i') {
            c++;
            te->oht = OHT_FIBM;
        } else {
            te->oht = OHT_FIEEE;
        }
        maxoff = THEAD_BYTES - 4;
        break;

    default: return;
    }

    c++;
    if(*c == 0) return;

    te->off = atoi(c);

    if(te->off < 1 || te->off > maxoff) {
        te->off = -1;
        return;
    }

    while(*c && isdigit(*c)) c++;
    while(*c && isspace(*c)) c++;
    if(*c != 0) {
        te->off = -1;
        return;
    }
    te->off -= 1;

    TRACE(("Using fast trace evaluation for [%s]: offset=%d", expr, te->off));
}

static double offset_only_evaluation(trc_eval* te, const byte* b)
{
    b += te->off;

    switch(te->oht) {
    case OHT_INT:
        return (double)((int)B2I32(b[0], b[1], b[2], b[3]));
    case OHT_UINT:
        return (double)((unsigned int)B2I32(b[0], b[1], b[2], b[3]));
    case OHT_SHORT:
        return (double)((short)B2I16(b[0], b[1]));
    case OHT_USHORT:
        return (double)((unsigned short)B2I16(b[0], b[1]));
    case OHT_BYTE:
        return (double)((char)b[0]);
    case OHT_UBYTE:
        return (double)((unsigned char)b[0]);
    case OHT_FIBM: {
        union {
            unsigned int i;
            float f;
        }i2f;
        unsigned int ival = (unsigned int)B2I32(b[0], b[1], b[2], b[3]);
        ibm_to_ieee(&ival, &i2f.i, 1, 1);
        return (double)i2f.f;
    }
    case OHT_FIEEE: {
        union {
            unsigned int val;
            float f;
        }i2f;
        i2f.val = (unsigned int)B2I32(b[0], b[1], b[2], b[3]);
        return (double)i2f.f;
    }

    default:
        ERROR(("Internal trace evaluation error"));
    }

    return 0.0;
}
