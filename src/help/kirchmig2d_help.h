#ifndef KIRCHMIG2D_HELP_H
#define KIRCHMIG2D_HELP_H

#include <iostream>
#include <string>

namespace kirchmig2d_help {

const std::string PROGRAM_TITLE = "2-D prestack Kirchhoff depth migration/modeling - kirchmig2d";
const std::string SEPARATOR_LINE = "========================================================";
const std::string DESCRIPTION = R"HELP(
kirchmig2d applies a 2-D prestack Kirchhoff operator. With adj=1, it reads
seismic data and writes a migrated depth image. With adj=0, it reads a depth
model/image and writes modeled prestack data.

The operator uses source and receiver traveltime tables plus their derivative
tables for interpolation and anti-aliasing. The image grid is obtained from
the source traveltime table stable=<file>.
)HELP";
const std::string USAGE =
    "kirchmig2d data=<data> mig=<image> stable=<stable> sderiv=<sderiv> rtable=<rtable> rderiv=<rderiv> [parameters]";
const std::string PARAMETERS = R"HELP(
Required parameters in all modes:
  data=<file>
      Seismic data file. If adj=1, this file is opened for reading and must
      contain input prestack data. If adj=0, this file is opened for writing
      and receives the modeled prestack data.

  mig=<file>
      Migrated image/model file. If adj=1, this file is opened for writing
      and receives the migrated image. If adj=0, this file is opened for
      reading and supplies the input model/image.

  stable=<file>
      Source traveltime table. This file defines the image grid. Required
      headers are n1, o1, d1 for depth; n2, o2, d2 for lateral position;
      and n3, o3, d3 for source-table coordinate. The code reads this table
      as n3 slices, each containing n1*n2 samples.

  sderiv=<file>
      Source traveltime derivative table. It must have the same dimensions
      and ordering as stable=<file>. The samples are used as stablex in the
      anti-aliasing estimate and in Hermite interpolation.

  rtable=<file>
      Receiver traveltime table. The code reads n3, o3, and d3 from this
      file as the receiver-table coordinate axis. The n1/n2 image-grid
      dimensions are assumed to be the same as stable=<file>; the program
      reads n3*n1*n2 samples.

  rderiv=<file>
      Receiver traveltime derivative table. It must have the same dimensions
      and ordering as rtable=<file>. The samples are used as rtablex in the
      anti-aliasing estimate and in Hermite interpolation.

Required headers when adj=1:
  data=<file> is read as input seismic data and must contain:
    n1, o1, d1  -> time axis: nt, t0, dt
    n2, o2, d2  -> offset or receiver axis: nh, h0, dh
    n3, o3, d3  -> shot axis: ns, s0, ds

Additional required parameters when adj=0:
  nt=<int>
      Number of time samples in the modeled data.

  dt=<float>
      Time sampling interval of the modeled data.

  dh=<float>
      Sampling interval of axis 2 in the modeled data. The current code
      requires dh even when nh defaults to 1; if nh=1, dh is reset to 0.0
      after parsing.

  ds=<float>
      Shot sampling interval of the modeled data. The current code requires
      ds even when ns defaults to 1; if ns=1, ds is reset to 0.0 after
      parsing.

Optional parameters:
  adj=<0|1>  [default: 1]
      Adjoint/modeling switch. adj=1 performs migration: data is input and
      mig is output. adj=0 performs modeling: mig is input and data is
      output.

  cmp=<0|1>  [default: 0]
      Acquisition-geometry switch for axis 2. If cmp=1, axis 2 is treated
      as offset and the receiver coordinate is computed as receiver=s+h. If
      cmp=0, axis 2 is treated as receiver coordinate directly. This flag
      also controls the output axis label: "Offset" for cmp=1 and
      "Receiver" for cmp=0.

  nh=<int>  [default for adj=0: 1]
      Number of samples on data axis 2 in modeling mode. In migration mode,
      nh is read from data header n2 and this parameter is ignored.

  ns=<int>  [default for adj=0: 1]
      Number of shot samples in modeling mode. In migration mode, ns is read
      from data header n3 and this parameter is ignored.

  t0=<float>  [default for adj=0: 0.0]
      Time origin of modeled data. In migration mode, t0 is read from data
      header o1.

  h0=<float>  [default for adj=0: 0.0]
      Origin of data axis 2. If cmp=1, this is offset origin. If cmp=0, this
      is receiver-coordinate origin. In migration mode, h0 is read from data
      header o2.

  s0=<float>  [default for adj=0: 0.0]
      Shot-coordinate origin of modeled data. In migration mode, s0 is read
      from data header o3.

  tau=<float>  [default: 0.0]
      Static time shift in seconds. The imaging/modeling time is computed as
      ti = source_traveltime + receiver_traveltime + tau.

  aperture=<float>  [default: 90.0]
      Kirchhoff migration aperture in degrees. The code applies a cone-angle
      aperture test for both source and receiver branches at each image
      point. Smaller values reject steeper source/receiver rays.

  antialias=<float>  [default: 1.0]
      Anti-aliasing scale factor. The local anti-aliasing parameter passed to
      kirmig_pick is antialias*max(abs(stablex*ds), abs(rtablex*dh)). Larger
      values apply stronger anti-aliasing.

  cig=<0|1>  [default: 0]
      Common-image-gather output switch. If cig=1 in migration mode, the
      output image has n3=nh and stores one image gather for each input axis-2
      sample. If cig=0, all axis-2 contributions are stacked into a single
      image with n3=1. In modeling mode, the input mig file must be arranged
      consistently with this choice.

  type=<linear|partial|hermit>  [default: hermit]
      Interpolation type for source and receiver traveltime tables. Only the
      first character is inspected by the code: 'l' selects linear
      interpolation, 'p' selects partial interpolation, and 'h' selects
      Hermite interpolation. Hermite interpolation uses both traveltime and
      derivative tables.

Output headers:
  If adj=1, mig=<file> is written on the image grid defined by stable:
    axis 1: Depth,   n1=nz, o1=z0, d1=dz
    axis 2: Lateral, n2=nx, o2=x0, d2=dx
    axis 3: n3=nh when cig=1, otherwise n3=1

  If adj=0, data=<file> is written with:
    axis 1: Time, n1=nt, o1=t0, d1=dt
    axis 2: Offset or Receiver, n2=nh, o2=h0, d2=dh
    axis 3: Shot, n3=ns, o3=s0, d3=ds

Processing notes:
  In migration mode, each input trace is double-integrated with doubint(...)
  before Kirchhoff summation. In modeling mode, the modeled trace is also
  double-integrated before being written. OpenMP is used when the code is
  compiled with SE_USE_OMP.

Boolean parameters:
  The current source code reads logical options with se_get_par_int(...).
  Therefore 0/1 values are the safest form on the command line.
)HELP";
const std::string EXAMPLE = R"HELP(
Examples:
  # Adjoint mode: prestack Kirchhoff depth migration
  kirchmig2d data=data.rsf mig=image.rsf \
             stable=stable.rsf sderiv=sderiv.rsf \
             rtable=rtable.rsf rderiv=rderiv.rsf \
             adj=1 cmp=0 aperture=90 antialias=1.0 cig=0 type=hermit

  # Modeling mode: generate prestack data from an input image/model
  kirchmig2d data=modeled.rsf mig=image.rsf \
             stable=stable.rsf sderiv=sderiv.rsf \
             rtable=rtable.rsf rderiv=rderiv.rsf \
             adj=0 cmp=0 nt=1500 dt=0.004 nh=200 dh=10 ns=80 ds=25 \
             t0=0 h0=0 s0=0 aperture=90 antialias=1.0 type=hermit
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

} // namespace kirchmig2d_help

#endif
