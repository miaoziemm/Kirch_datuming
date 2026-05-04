#ifndef EIKODS2D_HELP_H
#define EIKODS2D_HELP_H

#include <iostream>
#include <string>

namespace eikods2d_help {

const std::string PROGRAM_TITLE = "Eikonal traveltime and derivative solver - eikods2d";
const std::string SEPARATOR_LINE = "========================================================";
const std::string DESCRIPTION = R"HELP(
eikods2d computes first-arrival traveltime fields from an input velocity
or slowness-squared model and also writes traveltime-derivative auxiliary
volumes required by the Kirchhoff migration/datuming workflow.

The input model is read from in=<file>. Axis 1 and axis 2 are required;
axis 3 is optional and defaults to one sample. The code internally supports
an n1 x n2 x n3 grid, although a 2-D model is represented by n3=1.
)HELP";
const std::string USAGE =
    "eikods2d in=<model> out=<time> tdl1=<dl1> tds1=<ds1> [parameters]";
const std::string PARAMETERS = R"HELP(
Required parameters:
  in=<file>
      Input velocity or slowness-squared model. The file must contain n1,
      n2, d1, and d2 headers. Header n3 is optional and defaults to 1;
      d3 defaults to d2; o1, o2, and o3 default to 0.0 if absent.

  out=<file>
      Output traveltime file. One traveltime volume is written for each
      source position. The spatial headers are copied from the input model.

  tdl1=<file>
      Output file for the first required traveltime-derivative volume dl1.
      It has the same spatial sampling as out and is written once per shot.

  tds1=<file>
      Output file for the first required traveltime-derivative volume ds1.
      It has the same spatial sampling as out and is written once per shot.

Optional parameters:
  vel=<0|1>  [default: 1]
      Input-model type flag. If vel=1, in=<file> is interpreted as velocity
      and is converted internally to slowness squared, 1/v^2. If vel=0, the
      input samples are assumed to already be slowness squared.

  order=<int>  [default: 2]
      Accuracy order passed to the eikonal solver. The implementation is
      intended for order=1 or order=2. The code does not perform additional
      range checking, so values other than 1 or 2 should be avoided.

  sweep=<0|1>  [default: 0]
      Solver-selection flag parsed by the program. In the current source
      code, sweep=0 calls the fast-marching eikods(...) solver. The sweep=1
      branch is not implemented in the main loop and therefore should not be
      used for production calculations.

  br1=<float>  [default: d1]
  br2=<float>  [default: d2]
  br3=<float>  [default: d3]
      Physical half-size of the constant-velocity box around the source for
      axes 1, 2, and 3. These values are converted to sample counts through
      b1=round(br1/d1), b2=round(br2/d2), and b3=round(br3/d3), unless b1,
      b2, or b3 is supplied explicitly. The resulting b* values are forced
      to be at least one sample.

  b1=<int>
  b2=<int>
  b3=<int>
      Constant-velocity box half-size around the source in samples along
      axes 1, 2, and 3. These parameters override br1, br2, and br3. If a
      corresponding plane-wave flag is enabled and b* is not supplied, the
      default b* for that axis becomes the full axis length.

  plane1=<0|1>  [default: 0]
  plane2=<0|1>  [default: 0]
  plane3=<0|1>  [default: 0]
      Plane-wave source flags for axes 1, 2, and 3. When planei=1 and bi is
      not explicitly specified, the constant-velocity source box on that
      axis is expanded to the full model size. These flags are passed to the
      eikonal solver as the plane-wave source configuration.

  shotfile=<file>
      Optional source-location file. If supplied, the file must have n1=3
      and n2 equal to the number of shots. For each shot, the three samples
      are read as z, y, and x coordinates, respectively. With shotfile, the
      single-shot parameters zshot, yshot, and xshot are ignored.

  zshot=<float>  [default: 0.0]
  yshot=<float>  [default: o2 + 0.5*(n2-1)*d2]
  xshot=<float>  [default: o3 + 0.5*(n3-1)*d3]
      Source coordinates used only when shotfile is not supplied. The code
      stores source coordinates in the order z, y, x and passes them to the
      solver together with the corresponding model origins and spacings.

  tdl2=<file>
      Optional output file for the second traveltime-derivative volume dl2.
      If either tdl2 or tds2 is present, memory for second-derivative output
      is allocated.

  tds2=<file>
      Optional output file for the second traveltime-derivative volume ds2.
      If either tdl2 or tds2 is present, memory for second-derivative output
      is allocated.

  l=<int>  [default: 1]
      Source-perturbation direction index passed directly to eikods(...).
      Use the same convention as the derivative implementation in the
      SERECKIRCH library. The main program does not validate this value.

Header requirements for in=<file>:
  Required: n1, n2, d1, d2.
  Optional/defaulted: n3=1, d3=d2, o1=0, o2=0, o3=0.

Boolean parameters:
  The current source code reads logical options with se_get_par_int(...).
  Therefore 0/1 values are the safest form on the command line.
)HELP";
const std::string EXAMPLE = R"HELP(
Examples:
  # Single source at the default model center in axis 2 and axis 3
  eikods2d in=vel.rsf out=time.rsf tdl1=tdl1.rsf tds1=tds1.rsf \
           vel=1 order=2 br1=4 br2=4 br3=4 zshot=0

  # Multiple sources from a shot-location file; also write second derivatives
  eikods2d in=vel.rsf out=time.rsf shotfile=shots.rsf \
           tdl1=tdl1.rsf tds1=tds1.rsf tdl2=tdl2.rsf tds2=tds2.rsf \
           vel=1 order=2 l=1
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

} // namespace eikods2d_help

#endif
