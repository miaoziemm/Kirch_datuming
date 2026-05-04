#include "../include/se_pick_file.h"


const char* const PICKKEY_NHZ     = "nhz";
const char* const PICKKEY_OILINE  = "o_iline";
const char* const PICKKEY_OXLINE  = "o_xline";
const char* const PICKKEY_OZ      = "o_z";
const char* const PICKKEY_O4      = "o4";
const char* const PICKKEY_AZIMUTH = "azimuth";
const char* const PICKKEY_COLS    = "columns";
const char* const PICKKEY_AXES    = "axes";

const char* const PICKHKEY_LOAD     = "load";
const char* const PICKHKEY_COLOR    = "color";
const char* const PICKHKEY_FREEZE   = "freeze";
const char* const PICKHKEY_EXTENDED = "extended";

static const char* const NAME_HORIZON = "horizon";
static const char* const NAME_SET     = "set";

#define __PICK_LINE_MAXLEN (64*1024)

typedef struct {
    char* filename;
    FILE* f;
    int lineNo;
    hrz_surface_t* currHorizon;
    hrz_line_t* currSet;
    array currCoords;
    int inHorizon, inSet;
    double oi, ox, oz, o4;
    char** error;
    array* warns;
    hrz_model_t* model;
} __parser_state_t;

static void __init_state(__parser_state_t* state, char** error, array* warns);
static void __clear_state(__parser_state_t* state);


static int parseTypeLine(const char* line, const char* tag, char** type, __parser_state_t* state);
static int parseAttrLine(const char* line, char** key, char** val, __parser_state_t* state);
static int parsePointLine(const char* line, double* p, __parser_state_t* state);

static void parseOriginAttr(const char* key, double val, __parser_state_t* state);
static void applyOrigins(double* p, __parser_state_t* state);

static int checkAttrs(__parser_state_t* state);
static int isHorizonValid(__parser_state_t* state);
static int isSetValid(__parser_state_t* state);

static int checkRequiredAttr(const se_hash attrs, const char* key, const char* suffix, __parser_state_t* state);
static int checkOptionalAttr(const se_hash attrs, const char* key, const char* suffix, __parser_state_t* state);

static void setError(__parser_state_t* state, char* error);
static void addWarn(__parser_state_t* state, char* warn);

static void finalizeSet(__parser_state_t* state)
{
    int ncoords = array_size(state->currCoords);
    int i;

    state->currSet->npoints = ncoords / state->model->ndim;
    state->currSet->coords = (double *)malloc(ncoords * sizeof(double));
    for (i = 0; i < ncoords; ++i) {
        double* d = (double*) array_get_at(state->currCoords, i);
        state->currSet->coords[i] = *d;
    }
    array_clear(state->currCoords, 1);
    array_add(state->currHorizon->lines, state->currSet);
    state->currSet = (hrz_line_t *)malloc(sizeof(hrz_line_t));
    memset(state->currSet, 0, sizeof(hrz_line_t));
}

void __init_state(__parser_state_t* state, char** error, array* warns)
{
    state->filename = NULL;
    state->f = NULL;
    state->lineNo = 0;
    state->currHorizon = (hrz_surface_t *)malloc(sizeof(hrz_surface_t));
    init_hrz_surface(state->currHorizon);
    state->currSet = (hrz_line_t *)malloc(sizeof(hrz_line_t));
    state->currCoords = create_array();
    state->inHorizon = 0;
    state->inSet = 0;
    state->oi = state->ox = state->oz = state->o4 = 0.0;
    state->error = error;
    state->warns = warns;
    if (state->warns != NULL) {
        *(state->warns) = create_array();
    }
    state->model = (hrz_model_t *)malloc(sizeof(hrz_model_t));
    init_hrz_model(state->model);
}

void __clear_state(__parser_state_t* state)
{
    if (state->f != NULL) {
        if (fclose(state->f) != 0) {
            int err = errno;
            WARN(("Cannot close stream for file '%s' (current dir is '%s'): "
                  "errno %d - %s",
                  state->filename, get_cwd(), err, strerror(err)));
        }
    }
    free(state->filename);
    clear_hrz_surface(state->currHorizon);
    free(state->currHorizon);
    clear_hrz_line(state->currSet);
    free(state->currSet);
    destroy_array(state->currCoords, 1);
    clear_hrz_model(state->model);
}

hrz_model_t* parse_pick_file_print(const char* filename) {
    int i;
    hrz_model_t* model;
    char *error;
    array warns;

    model=parse_pick_file(filename, &error, &warns);

    // print warnings
    for (i = 0; i < array_size(warns); ++i) {
        const char* w = (const char*) array_get_at(warns, i);
        WARN((w));
    }
    destroy_array(warns, 1);

    if (model == NULL) {
        ERROR(("Parsing .pick file %s: %s",filename,error));
        free(error);
    }

    return model;
}

hrz_model_t* parse_pick_file(const char* filename, char** error, array* warns)
{
    __parser_state_t state;
    char line[__PICK_LINE_MAXLEN];
    int i;

    __init_state(&state, error, warns);

    /* open the file */
    state.filename = strdup(filename);
    state.f = fopen(filename, "r");
    if (state.f == NULL) {
        int err = errno;
        setError(&state, asprintf("Cannot open pick file '%s' (current dir is '%s'): errno %d - %s",
                                      filename, get_cwd(), err, strerror(err)));
        __clear_state(&state);
        return NULL;
    }

    /* parse each line */
    while (1) {
        /* read next line */
        ++state.lineNo;
        if (fgets(line, __PICK_LINE_MAXLEN, state.f) == NULL) {
            if (! feof(state.f)) {
                int err = errno;
                setError(&state, asprintf("Cannot read line from pick file '%s:%d': errno %d - %s",
                                              filename, state.lineNo, err, strerror(err)));
                __clear_state(&state);
                return NULL;
            }
            break;
        }

        strip_comment(line, '#');
        str_trim(line);
        if (strlen(line) == 0) {
            continue;
        }

        char* name;
        char* key;
        double p[4];

        name = NULL;
        if (parseTypeLine(line, NAME_HORIZON, &name, &state)) {
            if (! state.inHorizon) {
                /* validate the attribute */
                if (! checkAttrs(&state)) {
                    __clear_state(&state);
                    return NULL;
                }

                state.inHorizon = 1;
            } else {
                /* validate the last set */
                if (state.inSet) {
                    if (! isSetValid(&state)) {
                        __clear_state(&state);
                        free(name);
                        return NULL;
                    }
                    finalizeSet(&state);
                }

                /* validate the last horizon */
                if (! isHorizonValid(&state)) {
                    __clear_state(&state);
                    free(name);
                    return NULL;
                }
                array_add(state.model->surfaces, state.currHorizon);
                state.currHorizon = (hrz_surface_t *)malloc(sizeof(hrz_surface_t));
                init_hrz_surface(state.currHorizon);
            }
            free(state.currHorizon->name);
            state.currHorizon->name = name;
            state.inSet = 0;
            continue;
        }

        if (parseTypeLine(line, NAME_SET, &name, &state)) {
            /* validate the last set */
            if (state.inSet) {
                if (! isSetValid(&state)) {
                    __clear_state(&state);
                    free(name);
                    return NULL;
                }
                finalizeSet(&state);
            }

            state.inSet = 1;
            free(state.currSet->name);
            state.currSet->name = name;
            continue;
        }

        key = NULL;
        char *str_val = NULL;
        if (parseAttrLine(line, &key, &str_val, &state)) {
            if (! state.inHorizon) {
                ht_put(state.model->attrs, key, str_val);

                double dval;
                if (parse_double(str_val, &dval) == 0) {
                    parseOriginAttr(key, dval, &state);
                }
            } else if (! state.inSet) {
                ht_put(state.currHorizon->attrs, key, str_val);
            } else {
                addWarn(&state, asprintf("Line %d: Expected pick coords for horizon \"%s\" set \"%s\"",
                                             state.lineNo, state.currHorizon->name, state.currSet->name));
                free(key);
                free(str_val);
            }
            continue;
        }

        if (parsePointLine(line, p, &state)) {
            if (state.inSet) {
                /* swap columns 3 and 4 if we have 4 cols */
                if (state.model->ndim == 4) {
                    double tmp = p[2];
                    p[2] = p[3];
                    p[3] = tmp;
                }
                applyOrigins(p, &state);
                for (i = 0; i < state.model->ndim; ++i) {
                    double* d = (double *)malloc(sizeof(double));
                    *d = p[i];
                    array_add(state.currCoords, d);
                }
            } else {
                addWarn(&state, asprintf("Line %d: Unexpected pick coords outside horizon set", state.lineNo));
            }
            continue;
        }

        addWarn(&state, asprintf("Line %d: Unrecognized format", state.lineNo));
    }

    /* validate the last set */
    if (state.inSet && isSetValid(&state)) {
        finalizeSet(&state);
    }
    /* validate the last horizon */
    if (state.inHorizon && isHorizonValid(&state)) {
        array_add(state.model->surfaces, state.currHorizon);
        state.currHorizon = NULL;
    }

    /* validate horizons */
    int nHrz = array_size(state.model->surfaces);
    if (nHrz == 0) {
        addWarn(&state, strdup("No horizon defined"));

        /* validate the attributes */
        if (! checkAttrs(&state)) {
            __clear_state(&state);
            return NULL;
        }
    }
    int attrNhz = hrz_get_attr_i(state.model->attrs, PICKKEY_NHZ);
    if (attrNhz != nHrz) {
        addWarn(&state, asprintf("Expected %d horizon(s) (\"%s=%d\"), found %d",
                                     attrNhz, PICKKEY_NHZ, attrNhz, nHrz));
    }

    hrz_model_t* model = state.model;
    free(model->name);
    model->name = strdup(filename);

    /* cleanup */
    state.model = NULL;
    __clear_state(&state);

    /* success */
    return model;
}

int parseTypeLine(const char* line, const char* tag, char** type, __parser_state_t* state)
{
    char* s;

    /* split by first ' ' */
    int idx = str_indexof(line, ' ');
    if (idx < 0) {
        return 0;
    }

    char* tag_lower = strdup(tag);
    for (s = tag_lower; *s; ++s) {
        *s = (char)tolower(*s);
    }

    char* line_lower = strdup(line);
    for (s = line_lower; *s; ++s) {
        *s = (char)tolower(*s);
    }

    /* check tag (case insensitive) */
    int starts_with = str_startswith(line_lower, tag_lower);
    free(tag_lower);
    free(line_lower);
    if (! starts_with) {
        return 0;
    }

    /* get name (remove quotes if there and trim) */
    *type = strdup(line + idx + 1);
    int len = strlen(*type);
    if (len >= 2 && str_startswith(*type, "\"") && str_endswith(*type, "\"")) {
        memmove(*type, *type + 1, len - 2);
        (*type)[len - 2] = '\0';
    }
    str_trim(*type);
    if (strlen(*type) == 0) {
        addWarn(state, asprintf("Line %d: Empty %s name", state->lineNo, tag));
        free(*type);
        return 0;
    }

    return 1;
}

int parseAttrLine(const char* line, char** key, char** val, __parser_state_t* state)
{
    (void)state;
    /* split by '=' */
    int idx = str_indexof(line, '=');
    if (idx < 0) {
        return 0;
    }

    /* get key */
    *key = strndup(line, idx);

    /* get value */
    *val = strdup(line + idx + 1);

    return 1;
}

int parsePointLine(const char* line, double* p, __parser_state_t* state)
{
    int i;
    array sl = split2(line, ' ', SPLIT_NO_EMPTYFIELDS | SPLIT_TRIMFIELDS | SPLIT_USE_ISSPACE);
    if (array_size(sl) != state->model->ndim) {
        destroy_array(sl, 1);
        return 0;
    }

    memset(p, 0, 4 * sizeof(double));

    for (i = 0; i < array_size(sl); ++i) {
        const char* s = (const char*) array_get_at(sl, i);

        double dval;
        if (parse_double(s, &dval) == 0) {
            p[i] = dval;
        } else {
            destroy_array(sl, 1);
            return 0;
        }
    }

    destroy_array(sl, 1);
    return 1;
}

void parseOriginAttr(const char* key, double val, __parser_state_t* state)
{
    if (strcmp(key, PICKKEY_OILINE) == 0) {
        state->oi = val;
    } else if (strcmp(key, PICKKEY_OXLINE) == 0) {
        state->ox = val;
    } else if (strcmp(key, PICKKEY_OZ) == 0) {
        state->oz = val;
    } else if (strcmp(key, PICKKEY_O4) == 0) {
        state->o4 = val;
    }
}

void applyOrigins(double* p, __parser_state_t* state)
{
    if (state->model->ndim == 3) {
        p[0] += state->oi;
        p[1] += state->ox;
        p[2] += state->oz;
    } else if (state->model->ndim == 4) {
        p[0] += state->oi;
        p[1] += state->ox;
        p[2] += state->o4;
        p[3] += state->oz;
    }
}

int checkAttrs(__parser_state_t* state)
{
    if (! (checkRequiredAttr(state->model->attrs, PICKKEY_OILINE, "", state) &&
           checkRequiredAttr(state->model->attrs, PICKKEY_OXLINE, "", state) &&
           checkRequiredAttr(state->model->attrs, PICKKEY_OZ, "", state) &&
           checkRequiredAttr(state->model->attrs, PICKKEY_COLS, "", state))) {
        return 0;
    }
    checkOptionalAttr(state->model->attrs, PICKKEY_NHZ, "", state);
    checkOptionalAttr(state->model->attrs, PICKKEY_AZIMUTH, "", state);

    const char* sndim = hrz_get_attr(state->model->attrs, PICKKEY_COLS);
    if (parse_int(sndim, &state->model->ndim) != 0) {
        setError(state, asprintf("Invalid value for attribute '%s=%s': expected value >= 3", PICKKEY_COLS, sndim));
        return 0;
    }
    if (state->model->ndim < 3 || state->model->ndim > 4) {
        setError(state, asprintf("Invalid value for attribute '%s=%d': expected 3 or 4", PICKKEY_COLS, state->model->ndim));
        return 0;
    }

    return 1;
}

int isHorizonValid(__parser_state_t* state)
{
    char* suffix = asprintf("for horizon \"%s\"", state->currHorizon->name);
    checkOptionalAttr(state->currHorizon->attrs, PICKHKEY_LOAD, suffix, state);
    checkOptionalAttr(state->currHorizon->attrs, PICKHKEY_COLOR, suffix, state);
    checkOptionalAttr(state->currHorizon->attrs, PICKHKEY_FREEZE, suffix, state);
    checkOptionalAttr(state->currHorizon->attrs, PICKHKEY_EXTENDED, suffix, state);
    free(suffix);

    return 1;
}

int isSetValid(__parser_state_t* state)
{
    UNUSED(state);

    return 1;
}

int checkRequiredAttr(const se_hash attrs, const char* key, const char* suffix, __parser_state_t* state)
{
    if (! ht_contains_key(attrs, key)) {
        setError(state, asprintf("Missing required attribute '%s' %s", key, suffix));
        return 0;
    }
    return 1;
}

int checkOptionalAttr(const se_hash attrs, const char* key, const char* suffix, __parser_state_t* state)
{
    if (! ht_contains_key(attrs, key)) {
        addWarn(state, asprintf("Missing attribute '%s' %s", key, suffix));
        return 0;
    }
    return 1;
}

void setError(__parser_state_t* state, char* error)
{
    *(state->error) = error;
}

void addWarn(__parser_state_t* state, char* warn)
{
    if (state->warns != NULL) {
        array_add(*(state->warns), warn);
    } else {
        free(warn);
    }
}


int write_pick_file(const hrz_model_t* model, const char* filename, char** error)
{
    UNUSED(model);
    UNUSED(filename);
    UNUSED(error);

    // TODO

    return CODE_ERROR;
}
