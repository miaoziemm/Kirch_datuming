#!/usr/bin/env bash
set -euo pipefail

# Five-layer recursive Kirchhoff migration and ADCIG generation for the salt
# model. Put this script beside vsalt.rsf, data_salt_old.rsf, and
# par_salt_layer5.txt.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN_DIR="${BIN_DIR:-${SCRIPT_DIR}/build/bin}"
cd "${SCRIPT_DIR}"

nlayers=5
model_file="vsalt.rsf"
data_file="data_salt_old.rsf"
par_file="par_salt_layer5.txt"
layer_model_base="marmlayer.rsf"

nx=768
dx=0.01
cmp=0
tau=0.04
antialias=1
aperkir=40

# Half-opening-angle axis.
angle_min=0
angle_max=60
angle_step=1
angle_n=61
angle_interp=1
angle_taper_start=50
angle_taper_end=60
angle_smooth_rect=2
cig_x_index="${1:-384}"

# Source-centred V mute, one value per recursive layer. The boundary is
# t = intercept + slope*abs(receiver-source), and every sample with a smaller
# time is removed before a cosine recovery. Lateral coordinates use km, so the
# first-layer slope 50 s/km gives t=0.5 s at +/-0.01 km (+/-10 m).
# Decreasing slopes make the boundary progressively flatter with depth.
direct_mute_enable=(1 1 1 1 1)
direct_mute_invert=(0 0 0 0 0)
direct_mute_slope=(50.0 40.0 30.0 20.0 10.0)
direct_mute_intercept=(0.00 0.00 0.00 0.00 0.00)
direct_mute_taper=(0.05 0.06 0.07 0.08 0.10)
# Zero-based shot written before and after mute for per-layer quality control.
direct_mute_qc_shot=384

output_dir="./result/salt"
stack_output="${output_dir}/mig_salt_reckir_nlayer5_ADCIG_stack.rsf"
layer_scale_stack="1.0,1.5,1.4,2.4,2.8"
layer_scale_adcig="1.0,1.0,1.0,1.0,1.0"
stitch_blend=15
scalex_scale=1.3
scalex_aper=215

combined_adcig="${output_dir}/mig_salt_reckir_nlayer5_ADCIG_combined.rsf"
clean_adcig="${output_dir}/mig_salt_reckir_nlayer5_ADCIG_0_60.rsf"
single_raw="${output_dir}/mig_salt_reckir_nlayer5_ADCIG_x${cig_x_index}_0_60_raw.rsf"
single_smooth="${output_dir}/mig_salt_reckir_nlayer5_ADCIG_x${cig_x_index}_0_60_smooth.rsf"

fail() { echo "ERROR: $*" >&2; exit 1; }
require_file() { [[ -f "$1" ]] || fail "Missing input file: $1"; }
require_executable() { [[ -x "$1" ]] || fail "Missing executable: $1"; }
require_command() { command -v "$1" >/dev/null 2>&1 || fail "Missing Madagascar command: $1"; }
rsf_get_number()
{
    local file="$1" key="$2" value
    value="$(grep -Eo "(^|[[:space:]])${key}[[:space:]]*=[[:space:]]*[-+0-9.eE]+" "${file}" \
        | head -n 1 | sed -E "s/.*${key}[[:space:]]*=[[:space:]]*//")" || true
    printf '%s' "${value}"
}
remove_rsf()
{
    local file="$1"
    if command -v sfrm >/dev/null 2>&1 && [[ -f "${file}" ]]; then
        sfrm "${file}" >/dev/null 2>&1 || true
    fi
    rm -f "${file}" "${file}@"
}

require_file "${model_file}"
require_file "${data_file}"
require_file "${par_file}"
for executable in layervel shotgen eikods2d kirchmig2d_cig greensol_auto kirchdat2d_freq layercom layercom_cig extract_adcig; do
    require_executable "${BIN_DIR}/${executable}"
done
for command_name in sfstack sfin sfwindow sfsmooth; do
    require_command "${command_name}"
done

[[ "${cig_x_index}" =~ ^[0-9]+$ ]] || fail "cig_x_index must be a nonnegative integer"
(( cig_x_index < nx )) || fail "cig_x_index must be in [0,$((nx - 1))]"
[[ "${direct_mute_qc_shot}" =~ ^[0-9]+$ ]] || fail "direct_mute_qc_shot must be a nonnegative integer"
(( direct_mute_qc_shot < nx )) || fail "direct_mute_qc_shot must be in [0,$((nx - 1))]"
(( angle_n == (angle_max - angle_min) / angle_step + 1 )) || fail "Inconsistent angle axis"
for values in direct_mute_enable direct_mute_invert direct_mute_slope direct_mute_intercept direct_mute_taper; do
    declare -n array_ref="${values}"
    (( ${#array_ref[@]} == nlayers )) || fail "${values} must contain ${nlayers} values"
done
for ((i=1; i<nlayers; ++i)); do
    awk -v upper="${direct_mute_slope[i-1]}" -v lower="${direct_mute_slope[i]}" \
        'BEGIN {exit !(lower <= upper)}' || \
        fail "direct_mute_slope must not increase with depth (layers ${i} and $((i+1)))"
done
mkdir -p "${output_dir}"

build_tables()
{
    local layer="$1"
    "${BIN_DIR}/shotgen" shotfile="sht_${layer}.rsf" nshot="${nx}" oy=0 dy="${dx}"
    "${BIN_DIR}/shotgen" shotfile="rcv_${layer}.rsf" nshot="${nx}" oy=0 dy="${dx}"
    "${BIN_DIR}/eikods2d" \
        in="marmlayer_${layer}.rsf" out="time${layer}s.rsf" \
        shotfile="sht_${layer}.rsf" tdl1="tdl${layer}s.rsf" \
        tds1="tds${layer}s.rsf" b1=2 b2=2
    cp "time${layer}s.rsf" "time${layer}r.rsf"
    cp "tds${layer}s.rsf" "tds${layer}r.rsf"
    cp "tdl${layer}s.rsf" "tdl${layer}r.rsf"
}

make_layer_adcig()
{
    local layer="$1" input_data="$2" trace_aperture="$3"
    local index=$((layer - 1))

    printf '[Layer %d] mute: slope=%s s/km, intercept=%s s, taper=%s s\n' \
        "${layer}" "${direct_mute_slope[index]}" \
        "${direct_mute_intercept[index]}" "${direct_mute_taper[index]}"

    "${BIN_DIR}/kirchmig2d_cig" \
        seismic_data="${input_data}" migration="adcig_${layer}.rsf" \
        stable="time${layer}s.rsf" sderiv="tds${layer}s.rsf" \
        rtable="time${layer}r.rsf" rderiv="tds${layer}r.rsf" \
        adj=1 cig=1 mode=angle \
        angle_min="${angle_min}" angle_max="${angle_max}" \
        angle_step="${angle_step}" angle_n="${angle_n}" \
        angle_interp="${angle_interp}" \
        angle_taper_start="${angle_taper_start}" angle_taper_end="${angle_taper_end}" \
        direct_mute="${direct_mute_enable[index]}" \
        direct_mute_invert="${direct_mute_invert[index]}" \
        direct_mute_slope="${direct_mute_slope[index]}" \
        direct_mute_intercept="${direct_mute_intercept[index]}" \
        direct_mute_taper="${direct_mute_taper[index]}" \
        direct_mute_qc_shot="${direct_mute_qc_shot}" \
        direct_mute_qc_before="${output_dir}/mute_layer${layer}_shot${direct_mute_qc_shot}_before.rsf" \
        direct_mute_qc_after="${output_dir}/mute_layer${layer}_shot${direct_mute_qc_shot}_after.rsf" \
        cmp="${cmp}" tau="${tau}" antialias="${antialias}" \
        aperture="${aperkir}" aperture_trace="${trace_aperture}"

    sfstack axis=3 norm=n < "adcig_${layer}.rsf" > "mig_${layer}.rsf"
}

build_green_and_extrapolate()
{
    local layer="$1" input_data="$2" output_data="$3"
    "${BIN_DIR}/greensol_auto" ttabel_file="time${layer}s.rsf" \
        tgreen_file="time${layer}s_green.rsf" model_file="marmlayer_${layer}.rsf"
    "${BIN_DIR}/greensol_auto" ttabel_file="time${layer}r.rsf" \
        tgreen_file="time${layer}r_green.rsf" model_file="marmlayer_${layer}.rsf"
    "${BIN_DIR}/kirchdat2d_freq" \
        input_file="${input_data}" output_file="${output_data}" \
        aperture=200 data_taper=20 taper=0 length=0.05 \
        sgreen_file="time${layer}s_green.rsf" rgreen_file="time${layer}r_green.rsf" \
        model_file="marmlayer_${layer}.rsf" cmp="${cmp}" antialias=1 use_bf=1
}

"${BIN_DIR}/layervel" input_file="${model_file}" \
    output_file_base="${layer_model_base}" layer_txt="${par_file}"

for ((layer=1; layer<=nlayers; ++layer)); do
    echo "[Layer ${layer}/${nlayers}] Build traveltimes and ADCIG"
    build_tables "${layer}"
    if (( layer == 1 )); then input_data="${data_file}"; trace_aperture=200
    else input_data="rdata_$((layer - 1)).rsf"; trace_aperture=160
    fi
    make_layer_adcig "${layer}" "${input_data}" "${trace_aperture}"
    if (( layer < nlayers )); then
        build_green_and_extrapolate "${layer}" "${input_data}" "rdata_${layer}.rsf"
    fi
done

"${BIN_DIR}/layercom" \
    full_model="${model_file}" layer_model_base="${layer_model_base}" \
    layer_image_base="mig.rsf" layer_txt="${par_file}" \
    output_file="${stack_output}" scale="${layer_scale_stack}" \
    stitch_blend="${stitch_blend}" scalex_scale="${scalex_scale}" \
    scalex_aper="${scalex_aper}"

# Validate that all layer ADCIGs use the same angle axis before combining.
echo "[ADCIG 1/3] Validate and combine the five layer ADCIG volumes"
first_n3=""; first_o3=""; first_d3=""
for ((layer=1; layer<=nlayers; ++layer)); do
    adcig_file="adcig_${layer}.rsf"
    sfin "${adcig_file}" >/dev/null || fail "Cannot read ${adcig_file} or its binary data"
    n3="$(rsf_get_number "${adcig_file}" n3)"
    o3="$(rsf_get_number "${adcig_file}" o3)"; o3="${o3:-0}"
    d3="$(rsf_get_number "${adcig_file}" d3)"; d3="${d3:-1}"
    [[ -n "${n3}" ]] || fail "Cannot read n3 from ${adcig_file}"
    if [[ -z "${first_n3}" ]]; then
        first_n3="${n3}"; first_o3="${o3}"; first_d3="${d3}"
    else
        [[ "${n3}" == "${first_n3}" ]] || fail "Angle n3 mismatch in ${adcig_file}"
        awk -v a="${o3}" -v b="${first_o3}" \
            'BEGIN {d=a-b; if(d<0)d=-d; exit !(d<1e-6)}' || fail "Angle o3 mismatch in ${adcig_file}"
        awk -v a="${d3}" -v b="${first_d3}" \
            'BEGIN {d=a-b; if(d<0)d=-d; exit !(d<1e-6)}' || fail "Angle d3 mismatch in ${adcig_file}"
    fi
done

# Determine the available samples in the requested 0--60 degree interval.
read -r clean_f3 clean_n3 clean_o3 clean_d3 <<< "$(awk \
    -v n="${first_n3}" -v o="${first_o3}" -v d="${first_d3}" \
    -v amin="${angle_min}" -v amax="${angle_max}" '
    BEGIN {
      if(n<=0 || d<=0) exit 2;
      f=int((amin-o)/d+0.5); l=int((amax-o)/d+0.5);
      if(f<0)f=0; if(l>n-1)l=n-1; if(l<f)exit 3;
      printf "%d %d %.10g %.10g\n",f,l-f+1,o+f*d,d;
    }')" || fail "Cannot determine requested angle window"

for file in "${combined_adcig}" "${clean_adcig}" "${single_raw}" "${single_smooth}"; do
    remove_rsf "${file}"
done

"${BIN_DIR}/layercom_cig" \
    full_model="${model_file}" layer_image_base="adcig.rsf" \
    layer_txt="${par_file}" output_file="${combined_adcig}" \
    nlayers="${nlayers}" nangle="${first_n3}" scale="${layer_scale_adcig}" \
    stitch_blend="${stitch_blend}" scalex_scale="${scalex_scale}" \
    scalex_aper="${scalex_aper}"
sfin "${combined_adcig}" >/dev/null || fail "Combined ADCIG is unreadable"

echo "[ADCIG 2/3] Select the available 0--60 degree angle interval"
sfwindow f3="${clean_f3}" n3="${clean_n3}" \
    < "${combined_adcig}" > "${clean_adcig}"
sfin "${clean_adcig}" >/dev/null || fail "Windowed ADCIG is unreadable"

echo "[ADCIG 3/3] Extract ix=${cig_x_index} and smooth only along angle"
"${BIN_DIR}/extract_adcig" \
    input="${clean_adcig}" output="${single_raw}" ix="${cig_x_index}"
sfsmooth rect1=1 rect2="${angle_smooth_rect}" repeat=1 \
    < "${single_raw}" > "${single_smooth}"
for file in "${single_raw}" "${single_smooth}"; do
    [[ -s "${file}" && -s "${file}@" ]] || fail "Missing output: ${file}"
    sfin "${file}" >/dev/null || fail "Unreadable output: ${file}"
done

echo "Layered migration and ADCIG processing finished."
echo "Stacked image : ${stack_output}"
echo "Combined ADCIG: ${combined_adcig}"
echo "0--60 ADCIG   : ${clean_adcig}"
echo "Raw gather    : ${single_raw}"
echo "Smooth gather : ${single_smooth}"
