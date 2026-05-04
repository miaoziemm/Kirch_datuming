#ifndef SE_FS_SEGY_HEADER_H
#define SE_FS_SEGY_HEADER_H

#define EBC_COLS 80
#define EBC_ROWS 40
#define EBC_BYTES (EBC_COLS*EBC_ROWS) /* 3200 */
#define BHEAD_BYTES 400
#define THEAD_BYTES 240
#define NTRACES(fs, ns, bps) ( ( (fs) - EBC_BYTES - BHEAD_BYTES ) /     \
                               ( THEAD_BYTES + (ns) * (bps)) )

#define MAX_TAPE_RECORD 65535

#define THEAD_NKEYS 79
#define BHEAD_NKEYS 27

/** SEGY - trace identification header (240 bytes). */
typedef struct trchead_s {

    /** Byte 1: trace sequence number within line */
    int tracl;

    /** Byte 5: trace sequence number within reel */
    int tracr;

    /** Byte 9: field record number */
    int fldr;

    /** Byte 13: trace number within field record */
    int tracf;

    /** Byte 17: energy source point number */
    int ep;

    /** Byte 21: CDP ensemble number */
    int cdp;

    /** Byte 25: trace number within CDP ensemble */
    int cdpt;

    /** Byte 29: trace identification code:
       - -1 = Other
       - 0 = Unknown
       - 1 = Seismic data
       - 2 = Dead
       - 3 = Dummy
       - 4 = Time break
       - 5 = Uphole
       - 6 = Sweep
       - 7 = Timing
       - 8 = Waterbreak
       - 9 = Near-field gun signature
       - 10 = Far-field gun signature
       - 11 = Seismic pressure sensor
       - 12 = Multicomponent seismic sensor - Vertical component
       - 13 = Multicomponent seismic sensor - Cross-line component
       - 14 = Multicomponent seismic sensor - In-line component
       - 15 = Rotated multicomponent seismic sensor - Vertical component
       - 16 = Rotated multicomponent seismic sensor - Transverse component
       - 17 = Rotated multicomponent seismic sensor - Radial component
       - 18 = Vibrator reaction mass
       - 19 = Vibrator baseplate
       - 20 = Vibrator estimated ground force
       - 21 = Vibrator reference
       - 22 = Time-velocity pairs
       - 23 … N = optional use,  (maximum N = 32,767) */
    short trid;

    /** Byte 31: number of vertically summed traces (see vscode in bhed structure) */
    short nvs;

    /** Byte 33: number of horizontally summed traces (see vscode in bhed structure) */
    short nhs;

    /** Byte 35: data use:
        - 1 = production
        - 2 = test */
    short duse;

    /** Byte 37: distance from center of the source point to the
        center of the receiver group (negative if opposite to
        direction in which line is shot). */
    int offset;

    /** Byte 41: Receiver group elevation (all elevations above the
        Vertical datum are positive and below are negative). The scalar
        in bytes 69-70 applies to these values. The units are feet or
        meters as specified in Binary File Header bytes 55-56. */
    int gelev;

    /** Byte 45: Source elevation from sea level (above sea level is
        positive). The scalar in bytes 69-70 applies to these
        values. The units are feet or meters as specified in Binary
        File Header bytes 55-56. */
    int selev;

    /** Byte 49: source depth (positive). The scalar in bytes 69-70
        applies to these values. The units are feet or meters as
        specified in Binary File Header bytes 55-56.*/
    int sdepth;        

    /** Byte 53: datum elevation at receiver group. The scalar in
        bytes 69-70 applies to these values. The units are feet or
        meters as specified in Binary File Header bytes 55-56. */
    int gdel;
    
    /** Byte 57: datum elevation at source. The scalar in bytes 69-70
        applies to these values. The units are feet or meters as
        specified in Binary File Header bytes 55-56. */
    int sdel;

    /** Byte 61: water depth at source. The scalar in bytes 69-70
        applies to these values. The units are feet or meters as
        specified in Binary File Header bytes 55-56.*/
    int swdep;        

    /** Byte 65: water depth at receiver group. The scalar in bytes
        69-70 applies to these values. The units are feet or meters as
        specified in Binary File Header bytes 55-56. */
    int gwdep;

    /** Byte 69: Scalar to be applied to all elevations and depths
        specified in Trace Header bytes 41‑68 to give the real value.
        Scalar = 1, +10, +100, +1000, or +10,000.  If positive, scalar
        is used as a multiplier; if negative, scalar is used as a
        divisor.*/
    short scalel;

    /** Byte 71: Scalar to be applied to all coordinates specified in
        Trace Header bytes 73‑88 and to bytes Trace Header 181-188 to
        give the real value.  Scalar = 1, +10, +100, +1000, or
        +10,000.  If positive, scalar is used as a multiplier; if
        negative, scalar is used as divisor.*/
    short scalco;

    /** Byte 73: X source coordinate. If the coordinate units are in
        seconds of arc, decimal degrees or DMS, the X values represent
        longitude and the Y values latitude. A positive value
        designates east of Greenwich Meridian or north of the equator
        and a negative value designates south or west. */
    int  sx;

    /** Byte 77: Y source coordinate. If the coordinate units are in
        seconds of arc, decimal degrees or DMS, the X values represent
        longitude and the Y values latitude. A positive value
        designates east of Greenwich Meridian or north of the equator
        and a negative value designates south or west. */
    int  sy;

    /** Byte 81: X group coordinate. If the coordinate units are in
        seconds of arc, decimal degrees or DMS, the X values represent
        longitude and the Y values latitude. A positive value
        designates east of Greenwich Meridian or north of the equator
        and a negative value designates south or west. */
    int  gx;

    /** Byte 85: Y group coordinate. If the coordinate units are in
        seconds of arc, decimal degrees or DMS, the X values represent
        longitude and the Y values latitude. A positive value
        designates east of Greenwich Meridian or north of the equator
        and a negative value designates south or west. */
    int  gy;

    /** Byte 89: coordinate units code for previous four entries.
        - 1 = length (meters or feet)
        - 2 = seconds of arc
        - 3 = Decimal degrees
        - 4 = Degrees, minutes, seconds (DMS)

        \note To encode +/-DDDMMSS bytes 89-90 equal 
              +/-DDD*104 + MM*102 + SS with bytes 71-72 set to 1; 
              To encode +/- DDDMMSS.ss bytes 89-90 equal 
              +/-DDD*106 + MM*104 + SS*102 with bytes 71-72 set to -100. */
    short counit;        

    short wevel;        /* weathering velocity */

    short swevel;        /* subweathering velocity */

    short sut;        /* uphole time at source */

    short gut;        /* uphole time at receiver group */

    short sstat;        /* source static correction */

    short gstat;        /* group static correction */

    short tstat;        /* total static applied */

    short laga;        /* lag time A, time in ms between end of 240-
                          byte trace identification header and time
                          break, positive if time break occurs after
                          end of header, time break is defined as
                          the initiation pulse which maybe recorded
                          on an auxiliary trace or as otherwise
                          specified by the recording system */

    short lagb;        /* lag time B, time in ms between the time break
                          and the initiation time of the energy source,
                          may be positive or negative */

    short delrt;        /* delay recording time, time in ms between
                           initiation time of energy source and time
                           when recording of data samples begins
                           (for deep water work if recording does not
                           start at zero time) */

    short muts;        /* mute time--start */

    short mute;        /* mute time--end */

    unsigned short ns;        /* number of samples in this trace */

    unsigned short dt;        /* sample interval; in micro-seconds */

    short gain;        /* gain type of field instruments code:
                          1 = fixed
                          2 = binary
                          3 = floating point
                          4 ---- N = optional use */

    short igc;        /* instrument gain constant */

    short igi;        /* instrument early or initial gain */

    short corr;        /* correlated:
                          1 = no
                          2 = yes */

    short sfs;        /* sweep frequency at start */

    short sfe;        /* sweep frequency at end */

    short slen;        /* sweep length in ms */

    short styp;        /* sweep type code:
                          1 = linear
                          2 = cos-squared
                          3 = other */

    short stas;        /* sweep trace length at start in ms */

    short stae;        /* sweep trace length at end in ms */

    short tatyp;        /* taper type: 1=linear, 2=cos^2, 3=other */

    short afilf;        /* alias filter frequency if used */

    short afils;        /* alias filter slope */

    short nofilf;        /* notch filter frequency if used */

    short nofils;        /* notch filter slope */

    short lcf;        /* low cut frequency if used */

    short hcf;        /* high cut frequncy if used */

    short lcs;        /* low cut slope */

    short hcs;        /* high cut slope */

    short year;        /* year data recorded */

    short day;        /* day of year */

    short hour;        /* hour of day (24 hour clock) */

    short minute;        /* minute of hour */

    short sec;        /* second of minute */

    short timbas;        /* time basis code:
                            1 = local
                            2 = GMT
                            3 = other */

    short trwf;        /* trace weighting factor, defined as 1/2^N
                          volts for the least sigificant bit */

    short grnors;        /* geophone group number of roll switch
                            position one */

    short grnofr;        /* geophone group number of trace one within
                            original field record */

    short grnlof;        /* geophone group number of last trace within
                            original field record */

    short gaps;        /* gap size (total number of groups dropped) */

    short otrav;        /* overtravel taper code:
                           1 = down (or behind)
                           2 = up (or ahead) */

    /** Crossline number. */
    int xline;   /* byte 181 */

    /** Inline number. */
    int iline;   /* byte 185 */

    /** sample spacing between traces */
    int td2;     /* byte 189 */

    /** first trace location */
    int tf2;     /* byte 193 */

    /** negative of power used for dynamic range compression */
    int ungpow;  /* byte 197 */

    /** reciprocal of scaling factor to normalize range */
    int unscale; /* byte 201 */

    /** number of traces */
    int ntr;     /* byte 205 */

    /** mark selected traces */
    short mark;  /* byte 209 */

    short unass[15]; /* unassigned--NOTE: last entry causes 
                        a break in the word alignment, if we REALLY
                        want to maintain 240 bytes, the following
                        entry should be an odd number of short/UINT2
                        OR do the insertion above the "mark" keyword
                        entry */

} trchead;


typedef struct binhead_s {        /* bhed - binary header */

    int jobid;        /* job identification number */

    int lino;        /* line number (only one line per reel) */

    int reno;        /* reel number */

    short ntrpr;        /* number of data traces per record */

    short nart;        /* number of auxiliary traces per record */

    unsigned short hdt;  /* sample interval in micro secs for this reel */

    short dto;        /* same for original field recording */

    short hns;        /* number of samples per trace for this reel */

    short nso;        /* same for original field recording */

    short format;        /* data sample format code:
                            1 = floating point (4 bytes)
                            2 = fixed point (4 bytes)
                            3 = fixed point (2 bytes)
                            4 = fixed point w/gain code (4 bytes) */

    short fold;        /* CDP fold expected per CDP ensemble */

    short tsort;        /* trace sorting code: 
                           1 = as recorded (no sorting)
                           2 = CDP ensemble
                           3 = single fold continuous profile
                           4 = horizontally stacked */

    short vscode;        /* vertical sum code:
                            1 = no sum
                            2 = two sum ...
                            N = N sum (N = 32,767) */

    short hsfs;        /* sweep frequency at start */

    short hsfe;        /* sweep frequency at end */

    short hslen;        /* sweep length (ms) */

    short hstyp;        /* sweep type code:
                           1 = linear
                           2 = parabolic
                           3 = exponential
                           4 = other */

    short schn;        /* trace number of sweep channel */

    short hstas;        /* sweep trace taper length at start if
                           tapered (the taper starts at zero time
                           and is effective for this length) */

    short hstae;        /* sweep trace taper length at end (the ending
                           taper starts at sweep length minus the taper
                           length at end) */

    short htatyp;        /* sweep trace taper type code:
                            1 = linear
                            2 = cos-squared
                            3 = other */

    short hcorr;        /* correlated data traces code:
                           1 = no
                           2 = yes */

    short bgrcv;        /* binary gain recovered code:
                           1 = yes
                           2 = no */

    short rcvm;        /* amplitude recovery method code:
                          1 = none
                          2 = spherical divergence
                          3 = AGC
                          4 = other */

    short mfeet;        /* measurement system code:
                           1 = meters
                           2 = feet */

    short polyt;        /* impulse signal polarity code:
                           1 = increase in pressure or upward
                           geophone case movement gives
                           negative number on tape
                           2 = increase in pressure or upward
                           geophone case movement gives
                           positive number on tape */

    short vpol;        /* vibratory polarity code:
                          code        seismic signal lags pilot by
                          1        337.5 to  22.5 degrees
                          2         22.5 to  67.5 degrees
                          3         67.5 to 112.5 degrees
                          4        112.5 to 157.5 degrees
                          5        157.5 to 202.5 degrees
                          6        202.5 to 247.5 degrees
                          7        247.5 to 292.5 degrees
                          8        293.5 to 337.5 degrees */

    short hunass[170];        /* unassigned */

} binhead;

#define NOTHING()

#define SEGY_BHEADER_FIELDS_LIST(applyi, applys, applyus) NOTHING()     \
        applyi(jobid)     ,                                             \
        applyi(lino)      ,                                             \
        applyi(reno)      ,                                             \
        applys(ntrpr)     ,                                             \
        applys(nart)      ,                                             \
        applyus(hdt)      ,                                             \
        applys(dto)       ,                                             \
        applys(hns)       ,                                             \
        applys(nso)       ,                                             \
        applys(format)    ,                                             \
        applys(fold)      ,                                             \
        applys(tsort)     ,                                             \
        applys(vscode)    ,                                             \
        applys(hsfs)      ,                                             \
        applys(hsfe)      ,                                             \
        applys(hslen)     ,                                             \
        applys(hstyp)     ,                                             \
        applys(schn)      ,                                             \
        applys(hstas)     ,                                             \
        applys(hstae)     ,                                             \
        applys(htatyp)    ,                                             \
        applys(hcorr)     ,                                             \
        applys(bgrcv)     ,                                             \
        applys(rcvm)      ,                                             \
        applys(mfeet)     ,                                             \
        applys(polyt)     ,                                             \
        applys(vpol)


#define SEGY_HEADER_FIELDS_LIST1(applyi, applys, applyus) NOTHING()     \
        applyi(tracl)       ,                                           \
        applyi(tracr)       ,                                           \
        applyi(fldr)        ,                                           \
        applyi(tracf)       ,                                           \
        applyi(ep)          ,                                           \
        applyi(cdp)         ,                                           \
        applyi(cdpt)        ,                                           \
        applys(trid)        ,                                           \
        applys(nvs)         ,                                           \
        applys(nhs)         ,                                           \
        applys(duse)        ,                                           \
        applyi(offset)      ,                                           \
        applyi(gelev)       ,                                           \
        applyi(selev)       ,                                           \
        applyi(sdepth)      ,                                           \
        applyi(gdel)        ,                                           \
        applyi(sdel)        ,                                           \
        applyi(swdep)       ,                                           \
        applyi(gwdep)       ,                                           \
        applys(scalel)      ,                                           \
        applys(scalco)      ,                                           \
        applyi(sx)          ,                                           \
        applyi(sy)          ,                                           \
        applyi(gx)          ,                                           \
        applyi(gy)          ,                                           \
        applys(counit)      ,                                           \
        applys(wevel)       ,                                           \
        applys(swevel)      ,                                           \
        applys(sut)         ,                                           \
        applys(gut)         ,                                           \
        applys(sstat)       ,                                           \
        applys(gstat)       ,                                           \
        applys(tstat)       ,                                           \
        applys(laga)        ,                                           \
        applys(lagb)        ,                                           \
        applys(delrt)       ,                                           \
        applys(muts)        ,                                           \
        applys(mute)        ,                                           \
        applyus(ns)         ,                                           \
        applyus(dt)         ,                                           \
        applys(gain)        ,                                           \
        applys(igc)         ,                                           \
        applys(igi)         ,                                           \
        applys(corr)        ,                                           \
        applys(sfs)         ,                                           \
        applys(sfe)         ,                                           \
        applys(slen)        ,                                           \
        applys(styp)        ,                                           \
        applys(stas)        ,                                           \
        applys(stae)        ,                                           \
        applys(tatyp)       ,                                           \
        applys(afilf)       ,                                           \
        applys(afils)       ,                                           \
        applys(nofilf)      ,                                           \
        applys(nofils)      ,                                           \
        applys(lcf)         ,                                           \
        applys(hcf)         ,                                           \
        applys(lcs)         ,                                           \
        applys(hcs)         ,                                           \
        applys(year)        ,                                           \
        applys(day)         ,                                           \
        applys(hour)        ,                                           \
        applys(minute)      ,                                           \
        applys(sec)         ,                                           \
        applys(timbas)      ,                                           \
        applys(trwf)        ,                                           \
        applys(grnors)      ,                                           \
        applys(grnofr)      ,                                           \
        applys(grnlof)      ,                                           \
        applys(gaps)        ,                                           \
        applys(otrav)       ,                                           \
        applyi(xline)       ,                                           \
        applyi(iline)       ,                                           \
        applyi(td2)         ,                                           \
        applyi(tf2)         ,                                           \
        applyi(ungpow)      ,                                           \
        applyi(unscale)     ,                                           \
        applyi(ntr)         ,                                           \
        applys(mark)

#define SEGY_HEADER_FIELDS_LIST(apply) SEGY_HEADER_FIELDS_LIST1(apply, apply, apply)


#define SEGY_PTR_FIELD(f) *f

typedef struct segy_binhead_description_s {
    const char SEGY_BHEADER_FIELDS_LIST(SEGY_PTR_FIELD, SEGY_PTR_FIELD, SEGY_PTR_FIELD);
} segy_binhead_description_t;
extern const segy_binhead_description_t segy_binhead_description;

typedef struct segy_trchead_description_s {
    const char SEGY_HEADER_FIELDS_LIST1(SEGY_PTR_FIELD, SEGY_PTR_FIELD, SEGY_PTR_FIELD);
} segy_trchead_description_t;
extern const segy_trchead_description_t segy_trchead_description;



#endif