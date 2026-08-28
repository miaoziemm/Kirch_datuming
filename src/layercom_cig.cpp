#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace fs = std::filesystem;

struct RsfHeader {
    fs::path header_path;
    fs::path data_path;
    std::unordered_map<std::string, std::string> kv;

    int n1 = 0, n2 = 1, n3 = 1;
    double o1 = 0.0, o2 = 0.0, o3 = 0.0;
    double d1 = 1.0, d2 = 1.0, d3 = 1.0;
    int esize = 4;
    std::string data_format = "native_float";
    std::string label1 = "Depth";
    std::string label2 = "Lateral";
    std::string label3 = "Half-opening angle";
    std::string unit1;
    std::string unit2;
    std::string unit3 = "degree";
};

struct LayerRange {
    int start = 0;
    int end = 0;
    int should_datum_abs = -1;
};

static void fail(const std::string& msg) {
    throw std::runtime_error(msg);
}

static bool iequals(std::string a, std::string b) {
    std::transform(a.begin(), a.end(), a.begin(), [](unsigned char c){ return std::tolower(c); });
    std::transform(b.begin(), b.end(), b.begin(), [](unsigned char c){ return std::tolower(c); });
    return a == b;
}

static std::unordered_map<std::string, std::string>
parse_args(int argc, char** argv) {
    std::unordered_map<std::string, std::string> out;
    for (int i = 1; i < argc; ++i) {
        std::string s(argv[i]);
        const auto pos = s.find('=');
        if (pos == std::string::npos || pos == 0) {
            fail("Arguments must use key=value form; bad argument: " + s);
        }
        out[s.substr(0, pos)] = s.substr(pos + 1);
    }
    return out;
}

static std::string need_arg(const std::unordered_map<std::string, std::string>& a,
                            const std::string& key) {
    auto it = a.find(key);
    if (it == a.end() || it->second.empty()) fail("Need " + key + "=");
    return it->second;
}

static int get_int_arg(const std::unordered_map<std::string, std::string>& a,
                       const std::string& key, int def) {
    auto it = a.find(key);
    if (it == a.end()) return def;
    size_t used = 0;
    long long v = 0;
    try { v = std::stoll(it->second, &used); }
    catch (...) { fail("Invalid integer " + key + "=" + it->second); }
    if (used != it->second.size() || v < std::numeric_limits<int>::min() ||
        v > std::numeric_limits<int>::max()) {
        fail("Invalid integer " + key + "=" + it->second);
    }
    return static_cast<int>(v);
}

static double get_double_arg(const std::unordered_map<std::string, std::string>& a,
                             const std::string& key, double def) {
    auto it = a.find(key);
    if (it == a.end()) return def;
    size_t used = 0;
    double v = 0.0;
    try { v = std::stod(it->second, &used); }
    catch (...) { fail("Invalid number " + key + "=" + it->second); }
    if (used != it->second.size() || !std::isfinite(v)) {
        fail("Invalid number " + key + "=" + it->second);
    }
    return v;
}

static std::string strip_quotes(std::string s) {
    if (s.size() >= 2 && ((s.front() == '"' && s.back() == '"') ||
                          (s.front() == '\'' && s.back() == '\''))) {
        return s.substr(1, s.size() - 2);
    }
    return s;
}

static int kv_int(const std::unordered_map<std::string, std::string>& kv,
                  const std::string& key, int def, bool required = false) {
    auto it = kv.find(key);
    if (it == kv.end()) {
        if (required) fail("Missing RSF key " + key);
        return def;
    }
    size_t used = 0;
    int v = 0;
    try { v = std::stoi(it->second, &used); }
    catch (...) { fail("Invalid RSF integer " + key + "=" + it->second); }
    if (used != it->second.size()) fail("Invalid RSF integer " + key + "=" + it->second);
    return v;
}

static double kv_double(const std::unordered_map<std::string, std::string>& kv,
                        const std::string& key, double def) {
    auto it = kv.find(key);
    if (it == kv.end()) return def;
    size_t used = 0;
    double v = 0.0;
    try { v = std::stod(it->second, &used); }
    catch (...) { fail("Invalid RSF number " + key + "=" + it->second); }
    if (used != it->second.size() || !std::isfinite(v))
        fail("Invalid RSF number " + key + "=" + it->second);
    return v;
}

static std::string kv_string(const std::unordered_map<std::string, std::string>& kv,
                             const std::string& key, const std::string& def = "") {
    auto it = kv.find(key);
    return it == kv.end() ? def : it->second;
}

static RsfHeader read_rsf_header(const fs::path& header_path) {
    std::ifstream in(header_path);
    if (!in) fail("Cannot open RSF header: " + header_path.string());

    RsfHeader h;
    h.header_path = header_path;

    // Match key=value or key="quoted value" anywhere in a line.
    const std::regex re(R"(([A-Za-z_][A-Za-z0-9_]*)\s*=\s*("[^"]*"|[^\s]+))");
    std::string line;
    while (std::getline(in, line)) {
        for (std::sregex_iterator it(line.begin(), line.end(), re), end; it != end; ++it) {
            h.kv[(*it)[1].str()] = strip_quotes((*it)[2].str());
        }
    }

    h.n1 = kv_int(h.kv, "n1", 0, true);
    h.n2 = kv_int(h.kv, "n2", 1);
    h.n3 = kv_int(h.kv, "n3", 1);
    h.o1 = kv_double(h.kv, "o1", 0.0);
    h.o2 = kv_double(h.kv, "o2", 0.0);
    h.o3 = kv_double(h.kv, "o3", 0.0);
    h.d1 = kv_double(h.kv, "d1", 1.0);
    h.d2 = kv_double(h.kv, "d2", 1.0);
    h.d3 = kv_double(h.kv, "d3", 1.0);
    h.esize = kv_int(h.kv, "esize", 4);
    h.data_format = kv_string(h.kv, "data_format", "native_float");
    h.label1 = kv_string(h.kv, "label1", "Depth");
    h.label2 = kv_string(h.kv, "label2", "Lateral");
    h.label3 = kv_string(h.kv, "label3", "Half-opening angle");
    h.unit1 = kv_string(h.kv, "unit1", "");
    h.unit2 = kv_string(h.kv, "unit2", "");
    h.unit3 = kv_string(h.kv, "unit3", "degree");

    if (h.n1 <= 0 || h.n2 <= 0 || h.n3 <= 0)
        fail("Invalid RSF dimensions in " + header_path.string());
    if (h.esize != 4 || h.data_format != "native_float") {
        fail("Only data_format=\"native_float\", esize=4 is supported: " +
             header_path.string());
    }

    const std::string in_value = kv_string(h.kv, "in", "");
    if (in_value.empty()) fail("Missing in= in RSF header: " + header_path.string());

    fs::path p(in_value);
    if (p.is_absolute()) {
        h.data_path = p;
    } else {
        // Madagascar headers in this project use in="adcig_N.rsf@" next to the header.
        const fs::path beside_header = header_path.parent_path() / p;
        if (fs::exists(beside_header)) h.data_path = beside_header;
        else h.data_path = p;
    }
    return h;
}

static uint64_t sample_count(int n1, int n2, int n3) {
    if (n1 <= 0 || n2 <= 0 || n3 <= 0) fail("Invalid nonpositive dimension");
    return static_cast<uint64_t>(n1) * static_cast<uint64_t>(n2) *
           static_cast<uint64_t>(n3);
}

static std::vector<float> read_float_data(const RsfHeader& h) {
    const uint64_t count = sample_count(h.n1, h.n2, h.n3);
    const uint64_t bytes = count * sizeof(float);

    std::error_code ec;
    const auto actual = fs::file_size(h.data_path, ec);
    if (ec) fail("Cannot stat RSF binary: " + h.data_path.string());
    if (actual != bytes) {
        std::ostringstream oss;
        oss << "RSF binary size mismatch for " << h.header_path
            << ": expected " << bytes << " bytes, got " << actual
            << " bytes from " << h.data_path;
        fail(oss.str());
    }

    if (count > static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
        fail("RSF volume too large: " + h.header_path.string());
    std::vector<float> data(static_cast<size_t>(count));
    std::ifstream in(h.data_path, std::ios::binary);
    if (!in) fail("Cannot open RSF binary: " + h.data_path.string());
    in.read(reinterpret_cast<char*>(data.data()), static_cast<std::streamsize>(bytes));
    if (!in || static_cast<uint64_t>(in.gcount()) != bytes)
        fail("Short read from RSF binary: " + h.data_path.string());
    return data;
}

static void replace_all(std::string& s, const std::string& from, const std::string& to) {
    if (from.empty()) return;
    size_t pos = 0;
    while ((pos = s.find(from, pos)) != std::string::npos) {
        s.replace(pos, from.size(), to);
        pos += to.size();
    }
}

static std::vector<LayerRange> read_layer_ranges(const fs::path& file,
                                                  int full_n1) {
    std::ifstream in(file);
    if (!in) fail("Cannot open layer_txt=" + file.string());

    std::vector<LayerRange> layers;
    std::string line;
    int line_no = 0;
    while (std::getline(in, line)) {
        ++line_no;
        const auto hash = line.find('#');
        if (hash != std::string::npos) line.resize(hash);
        replace_all(line, "，", " ");
        for (char& c : line) if (c == ',' || c == ';' || c == '\t') c = ' ';

        std::istringstream iss(line);
        std::vector<std::string> tok;
        for (std::string t; iss >> t;) tok.push_back(t);
        if (tok.empty()) continue;
        if (tok.size() < 2) fail("layer_txt line " + std::to_string(line_no) + " needs start,end");

        LayerRange r;
        if (iequals(tok[0], "START")) r.start = 0;
        else {
            try { r.start = std::stoi(tok[0]); }
            catch (...) { fail("Invalid layer start at line " + std::to_string(line_no)); }
        }

        if (iequals(tok[1], "END")) {
            r.end = full_n1 - 1;
            r.should_datum_abs = -1;
        } else {
            try { r.end = std::stoi(tok[1]); }
            catch (...) { fail("Invalid layer end at line " + std::to_string(line_no)); }
            if (tok.size() >= 3 && !iequals(tok[2], "END")) {
                try { r.should_datum_abs = std::stoi(tok[2]); }
                catch (...) { fail("Invalid should_datum_abs at line " + std::to_string(line_no)); }
            }
        }
        layers.push_back(r);
    }

    if (layers.empty()) fail("No layers found in " + file.string());
    if (layers.front().start != 0) fail("First layer must start at depth sample 0");
    for (size_t i = 0; i < layers.size(); ++i) {
        const auto& r = layers[i];
        if (r.start < 0 || r.end < r.start || r.end >= full_n1)
            fail("Invalid layer range [" + std::to_string(r.start) + "," +
                 std::to_string(r.end) + "]");
        if (i > 0) {
            if (r.start > layers[i-1].end)
                fail("Layer " + std::to_string(i+1) + " does not overlap previous layer");
            if (r.end <= layers[i-1].end)
                fail("Layer " + std::to_string(i+1) + " must extend deeper than previous layer");
        }
    }
    if (layers.back().end != full_n1 - 1)
        fail("Last layer must end at full-model sample " + std::to_string(full_n1 - 1));
    return layers;
}

static std::vector<float> parse_scales(const std::string& s, size_t nlayers) {
    std::vector<float> out(nlayers, 1.0f);
    if (s.empty()) return out;
    std::string t = s;
    replace_all(t, "，", ",");
    for (char& c : t) if (c == ';') c = ',';
    std::vector<float> vals;
    std::stringstream ss(t);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        size_t used = 0;
        float v = 0.0f;
        try { v = std::stof(item, &used); }
        catch (...) { fail("Invalid scale value: " + item); }
        while (used < item.size() && std::isspace(static_cast<unsigned char>(item[used]))) ++used;
        if (used != item.size() || !std::isfinite(v)) fail("Invalid scale value: " + item);
        vals.push_back(v);
    }
    if (vals.size() == 1) std::fill(out.begin(), out.end(), vals[0]);
    else if (vals.size() == nlayers) out = vals;
    else fail("scale= must contain 1 or " + std::to_string(nlayers) + " values");
    return out;
}

static bool axis_match(double a, double b) {
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) <= 1.0e-6 * scale;
}

static fs::path build_layer_name(const fs::path& base, int layer) {
    const std::string ext = base.extension().string();
    const std::string stem = base.stem().string();
    const fs::path parent = base.parent_path();
    if (ext == ".rsf") return parent / (stem + "_" + std::to_string(layer) + ext);
    return parent / (base.filename().string() + "_" + std::to_string(layer));
}

static inline size_t idx3(int i1, int i2, int i3, int n1, int n2) {
    return (static_cast<size_t>(i3) * static_cast<size_t>(n2) +
            static_cast<size_t>(i2)) * static_cast<size_t>(n1) +
           static_cast<size_t>(i1);
}

static void write_rsf(const fs::path& output_header,
                      const RsfHeader& first,
                      const RsfHeader& model,
                      int n1, int n2, int n3,
                      const std::vector<float>& data) {
    if (!output_header.parent_path().empty())
        fs::create_directories(output_header.parent_path());

    const fs::path output_binary = fs::path(output_header.string() + "@");
    {
        std::ofstream out(output_binary, std::ios::binary | std::ios::trunc);
        if (!out) fail("Cannot create output binary: " + output_binary.string());
        const uint64_t bytes = static_cast<uint64_t>(data.size()) * sizeof(float);
        out.write(reinterpret_cast<const char*>(data.data()), static_cast<std::streamsize>(bytes));
        if (!out) fail("Failed writing output binary: " + output_binary.string());
    }

    const uint64_t expected_bytes = sample_count(n1, n2, n3) * sizeof(float);
    std::error_code ec;
    const uint64_t actual_bytes = fs::file_size(output_binary, ec);
    if (ec || actual_bytes != expected_bytes) {
        std::ostringstream oss;
        oss << "Output binary size check failed: expected " << expected_bytes
            << ", got " << (ec ? 0 : actual_bytes);
        fail(oss.str());
    }

    {
        std::ofstream h(output_header, std::ios::trunc);
        if (!h) fail("Cannot create output header: " + output_header.string());
        h << "in=\"" << output_binary.string() << "\"\n";
        h << "data_format=\"native_float\"\n";
        h << "esize=4\n";
        h << "n1=" << n1 << "\n";
        h << "n2=" << n2 << "\n";
        h << "n3=" << n3 << "\n";
        h << std::setprecision(12);
        h << "o1=" << model.o1 << "\n";
        h << "d1=" << model.d1 << "\n";
        h << "o2=" << model.o2 << "\n";
        h << "d2=" << model.d2 << "\n";
        h << "o3=" << first.o3 << "\n";
        h << "d3=" << first.d3 << "\n";
        h << "label1=\"" << first.label1 << "\"\n";
        h << "label2=\"" << first.label2 << "\"\n";
        h << "label3=\"" << first.label3 << "\"\n";
        if (!first.unit1.empty()) h << "unit1=\"" << first.unit1 << "\"\n";
        if (!first.unit2.empty()) h << "unit2=\"" << first.unit2 << "\"\n";
        h << "unit3=\"" << (first.unit3.empty() ? "degree" : first.unit3) << "\"\n";
        if (!h) fail("Failed writing output header: " + output_header.string());
    }

    // Re-read our own header and verify it points at exactly the file we wrote.
    const RsfHeader check = read_rsf_header(output_header);
    if (check.n1 != n1 || check.n2 != n2 || check.n3 != n3)
        fail("Output header dimension verification failed");
    const uint64_t check_bytes = sample_count(check.n1, check.n2, check.n3) * sizeof(float);
    const uint64_t check_actual = fs::file_size(check.data_path, ec);
    if (ec || check_actual != check_bytes)
        fail("Output RSF header/binary verification failed");
}

int main(int argc, char** argv) {
    try {
        const auto args = parse_args(argc, argv);
        const fs::path model_file = need_arg(args, "full_model");
        const fs::path layer_base = need_arg(args, "layer_image_base");
        const fs::path layer_txt = need_arg(args, "layer_txt");
        const fs::path output_file = need_arg(args, "output_file");

        const int expected_nlayers = get_int_arg(args, "nlayers", 0);
        const int expected_nangle = get_int_arg(args, "nangle", 0);
        const int stitch_blend = get_int_arg(args, "stitch_blend", 15);
        const int scalex_aper = get_int_arg(args, "scalex_aper", 0);
        const double scalex_scale = get_double_arg(args, "scalex_scale", 1.0);
        const std::string scale_arg = args.count("scale") ? args.at("scale") : "";

        if (stitch_blend < -1) fail("stitch_blend must be -1 or nonnegative");
        if (scalex_aper < 0) fail("scalex_aper must be nonnegative");
        if (!std::isfinite(scalex_scale)) fail("Invalid scalex_scale");

        const RsfHeader model = read_rsf_header(model_file);
        const int full_n1 = model.n1;
        const int full_n2 = model.n2;
        auto layers = read_layer_ranges(layer_txt, full_n1);
        if (expected_nlayers > 0 && expected_nlayers != static_cast<int>(layers.size()))
            fail("nlayers does not match layer_txt");

        const fs::path first_file = build_layer_name(layer_base, 1);
        const RsfHeader first = read_rsf_header(first_file);
        const int nangle = first.n3;
        if (expected_nangle > 0 && expected_nangle != nangle)
            fail("nangle does not match first ADCIG n3");
        if (first.n2 != full_n2)
            fail("Lateral size mismatch between full_model and first ADCIG");
        if (!axis_match(first.o2, model.o2) || !axis_match(first.d2, model.d2))
            fail("Lateral axis mismatch between full_model and first ADCIG");

        std::cout << "INFO: CIG output dimensions: n1=" << full_n1
                  << " n2=" << full_n2 << " n3=" << nangle << "\n";
        std::cout << "INFO: Angle axis: o3=" << first.o3
                  << " d3=" << first.d3 << " n3=" << nangle << "\n";

        const auto layer_scale = parse_scales(scale_arg, layers.size());
        std::vector<float> lateral_scale(full_n2, 1.0f);
        if (scalex_aper > 0) {
            for (int i2 = 0; i2 < full_n2; ++i2) {
                const int edge = std::min(i2, full_n2 - 1 - i2);
                if (edge < scalex_aper) {
                    const double taper = (scalex_aper == 1)
                        ? 1.0
                        : 1.0 - static_cast<double>(edge) / static_cast<double>(scalex_aper - 1);
                    lateral_scale[i2] = static_cast<float>(1.0 + (scalex_scale - 1.0) * taper);
                }
            }
        }

        const uint64_t out_count64 = sample_count(full_n1, full_n2, nangle);
        if (out_count64 > static_cast<uint64_t>(std::numeric_limits<size_t>::max()))
            fail("Output CIG too large");
        std::vector<float> output(static_cast<size_t>(out_count64), 0.0f);

        int previous_end = -1;
        for (size_t il = 0; il < layers.size(); ++il) {
            const fs::path image_file = build_layer_name(layer_base, static_cast<int>(il + 1));
            const RsfHeader image = read_rsf_header(image_file);
            const int start = layers[il].start;
            const int end = layers[il].end;
            const int expected_n1 = end - start + 1;

            if (image.n1 != expected_n1 || image.n2 != full_n2 || image.n3 != nangle) {
                std::ostringstream oss;
                oss << image_file << " dimensions are " << image.n1 << "," << image.n2 << ","
                    << image.n3 << "; expected " << expected_n1 << "," << full_n2 << "," << nangle;
                fail(oss.str());
            }
            if (!axis_match(image.o2, first.o2) || !axis_match(image.d2, first.d2) ||
                !axis_match(image.o3, first.o3) || !axis_match(image.d3, first.d3))
                fail("Lateral or angle axis mismatch in " + image_file.string());

            int overlap = 0;
            if (il > 0) {
                overlap = previous_end - start + 1;
                if (overlap <= 0) fail("No overlap for layer " + std::to_string(il + 1));
            }
            int blend_count = overlap;
            if (overlap > 0 && stitch_blend >= 0) {
                blend_count = std::min(stitch_blend, overlap);
            }
            const int keep_upper = overlap - blend_count;

            std::cout << "INFO: Layer " << (il + 1) << ": " << image_file.string()
                      << ", range=[" << start << "," << end << "]"
                      << ", overlap=" << overlap
                      << ", keep_upper=" << keep_upper
                      << ", blend=" << blend_count
                      << ", scale=" << layer_scale[il] << "\n";

            const auto layer_data = read_float_data(image);
            long long nonfinite = 0;

#ifdef _OPENMP
#pragma omp parallel for collapse(2) reduction(+:nonfinite)
#endif
            for (int i3 = 0; i3 < nangle; ++i3) {
                for (int i2 = 0; i2 < full_n2; ++i2) {
                    const float sx = (il > 0) ? lateral_scale[i2] : 1.0f;
                    const float scale = layer_scale[il] * sx;
                    for (int i1 = 0; i1 < image.n1; ++i1) {
                        const int g1 = start + i1;
                        float v = layer_data[idx3(i1, i2, i3, image.n1, full_n2)] * scale;
                        if (!std::isfinite(v)) { ++nonfinite; v = 0.0f; }
                        const size_t oi = idx3(g1, i2, i3, full_n1, full_n2);

                        if (overlap > 0 && g1 <= previous_end) {
                            const int k = g1 - start;
                            if (k < keep_upper || blend_count == 0) continue;
                            const int bi = k - keep_upper;
                            const float w = (blend_count > 1)
                                ? static_cast<float>(bi) / static_cast<float>(blend_count - 1)
                                : 1.0f;
                            output[oi] = output[oi] * (1.0f - w) + v * w;
                        } else {
                            output[oi] = v;
                        }
                    }
                }
            }
            if (nonfinite > 0)
                fail(image_file.string() + " contains non-finite samples after scaling");
            previous_end = end;
        }

        write_rsf(output_file, first, model, full_n1, full_n2, nangle, output);
        std::cout << "INFO: Finished CIG stitching: " << output_file.string() << "\n";
        std::cout << "INFO: Output binary: " << output_file.string() << "@\n";
        std::cout << "INFO: Output samples: " << out_count64
                  << ", bytes=" << out_count64 * sizeof(float) << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 1;
    }
}