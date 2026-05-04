#ifndef KIRCHDAT2D_HELP_H
#define KIRCHDAT2D_HELP_H

#include <iostream>
#include <string>

namespace kirchdat2d_help {

const std::string PROGRAM_TITLE = "2-D prestack Kirchhoff redatuming - kirchdat2d";
const std::string SEPARATOR_LINE = "========================================================";
const std::string DESCRIPTION = R"HELP(
kirchdat2d performs two-stage prestack Kirchhoff redatuming on input
common-shot gathers. The first stage applies receiver-side continuation
using rgreen=<file>. The second stage applies source-side continuation
using sgreen=<file>. The final result is written to out=<file>.

The input data are assumed to be organized as n1=time, n2=receiver/offset
index within each common-shot gather, and n3=shot index.
)HELP";
const std::string USAGE =
    "kirchdat2d in=<data> out=<data_out> sgreen=<sgreen> rgreen=<rgreen> sdatum=<zs> rdatum=<zr> [parameters]";
const std::string PARAMETERS = R"HELP(
Required parameters:
  in=<file>
      Input prestack data. The file must contain n1, n2, n3, o1, d1, o2,
      d2, o3, and d3 headers. In the code, n1=nt, n2=nh, n3=ns, o1=t0,
      d1=dt, o2=h0, d2=dh, o3=s0, and d3=ds.

  out=<file>
      Output redatumed prestack data. The output sample count is nt*nh*ns,
      i.e., the same logical dimensions as the input data.

  sgreen=<file>
      Source-side Green's-function or traveltime table. The program reads
      n1, o1, and d1 from this file as nsg, sg0, and dsg, allocates an
      nsg-by-nsg table, and reads nsg*nsg floating-point samples. The table
      is indexed by source-grid positions during the source-side continuation.

  rgreen=<file>
      Receiver-side Green's-function or traveltime table. The program reads
      n1, o1, and d1 from this file as nrg, rg0, and drg, allocates an
      nrg-by-nrg table, and reads nrg*nrg floating-point samples. The table
      is indexed by receiver-grid positions during the receiver-side
      continuation.

  sdatum=<float>
      Source datum depth or vertical continuation distance used in the
      source-side Kirchhoff weight. It appears in the factor
      sdatum*tau/(sdatum^2 + horizontal_distance^2), so it should be in the
      same length unit as ds, dh, sg0, and dsg.

  rdatum=<float>
      Receiver datum depth or vertical continuation distance used in the
      receiver-side Kirchhoff weight. It appears in the factor
      rdatum*tau/(rdatum^2 + horizontal_distance^2), so it should be in the
      same length unit as ds, dh, rg0, and drg.

Optional parameters:
  aperture=<int>  [default: 50]
      Kirchhoff summation aperture measured in trace/sample indices. In the
      receiver-side pass, traces from ih-aperture to ih+aperture are used.
      In the source-side pass, source indices within the corresponding
      aperture are used, with the internal jump factor applied when ds and
      dh differ.

  taper=<int>  [default: 10]
      Number of trace/sample indices used for linear tapering at both edges
      of the aperture. A larger value gives a smoother aperture edge but
      reduces the effective contribution near the aperture boundary.

  length=<float>  [default: 0.025]
      Temporal filter length in seconds passed to filt_init(dt,length). The
      time sampling dt is read from input header d1.

  interm=<file>
      Optional intermediate output after receiver-side continuation and
      before source-side continuation. This is useful for checking whether
      receiver redatuming alone is producing reasonable gathers.

  verb=<0|1>  [default: 0]
      Verbosity flag. If verb=1, the program reports progress for common-
      shot gathers and common-receiver gathers.

Header requirements for in=<file>:
  Required: n1, n2, n3, o1, d1, o2, d2, o3, d3.
  Axis convention used by the code:
    axis 1: time,     n1=nt, o1=t0, d1=dt
    axis 2: receiver or offset coordinate, n2=nh, o2=h0, d2=dh
    axis 3: shot coordinate, n3=ns, o3=s0, d3=ds

Green's-function table requirements:
  sgreen and rgreen are read as square dense tables. The code checks n1,
  o1, and d1 only, but the actual data file must contain n1*n1 samples for
  each table. Coordinate ranges must cover all source and receiver indices
  requested by the acquisition geometry; otherwise the program stops with
  "Source table too small" or "Receiver table too small".

Boolean parameters:
  The current source code reads logical options with se_get_par_int(...).
  Therefore 0/1 values are the safest form on the command line.
)HELP";
const std::string EXAMPLE = R"HELP(
Examples:
  kirchdat2d in=data_surface.rsf out=data_datum.rsf \
             sgreen=sgreen.rsf rgreen=rgreen.rsf \
             sdatum=200 rdatum=200 aperture=50 taper=10 length=0.025 verb=1

  # Save the receiver-side continuation result for diagnosis
  kirchdat2d in=data_surface.rsf out=data_datum.rsf interm=data_rec_only.rsf \
             sgreen=sgreen.rsf rgreen=rgreen.rsf sdatum=200 rdatum=200
)HELP";

inline void print_help() {
    std::cout << SEPARATOR_LINE << std::endl;
    std::cout << PROGRAM_TITLE << std::endl;
    std::cout << SEPARATOR_LINE << std::endl;
    std::cout << DESCRIPTION << std::endl << std::endl;
    std::cout << "Usage:" << std::endl;
    std::cout << "  " << USAGE << std::endl << std::endl;
    std::cout << PARAMETERS << std::endl;
    std::cout << EXAMPLE << std::endl;
    std::cout << SEPARATOR_LINE << std::endl;
}

} // namespace kirchdat2d_help

#endif
