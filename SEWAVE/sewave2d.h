#ifndef SEWAVE2D_H
#define SEWAVE2D_H

#include <string>
#include <vector>

namespace sewave {

struct Grid2D { int nz=0,nx=0; float dz=1,dx=1,z0=0,x0=0; std::vector<float> v; };
struct Data3D { int nt=0,nr=0,ns=0; float dt=0.001f,dr=1,ds=1,t0=0,r0=0,s0=0; std::vector<float> d; };
struct ForwardParams { float dt=0.001f; int nt=1000; float fdom=20; float sx=0, sz=0; float rz=0; int nr=0; float r0=0, dr=1; int ns=1; float ds=0; int nbc=100,L=30; float alpha=1.0f; int type_compute_laplace=1; int type=0; bool visco=false; float q=1000; float f0=20,Omega0=0,fc=120; int order=2; bool flag_smooth=true, flag_homo=false; int progress_interval=0; };
struct ImageParams { bool cmp=true; bool compensate=false; int type=0; bool visco=false; float q=1000; int nbc=100,L=30; float alpha=1.0f; int type_compute_laplace=1; float sz=0, rz=0; int shot_begin=0, shot_end=-1; int amp_compensation_sign=1; float fdom=20, f0=20,Omega0=0,fc=120; int order=2; bool flag_smooth=true, flag_homo=false; int progress_interval=0; int max_parallel_shots=1; bool in_memory_snapshots=false; int snapshot_interval=1; std::vector<int> debug_snapshot_steps; std::string debug_wavefield_prefix; };

Grid2D read_rsf2d(const std::string& path);
Data3D read_rsf3d(const std::string& path);
void write_rsf2d(const std::string& path, const Grid2D& g, const std::vector<float>& data, const char* label="Amplitude");
void write_rsf3d(const std::string& path, const Data3D& d);
Data3D forward(const Grid2D& vel, const ForwardParams& p, const Grid2D* qmodel=nullptr);
std::vector<float> rtm_image(const Grid2D& vel, const Data3D& data, const ImageParams& p, const Grid2D* qmodel=nullptr);

} // namespace sewave
#endif
