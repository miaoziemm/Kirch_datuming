#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>

namespace fs = std::filesystem;

static std::string trim(const std::string& s) {
    const auto b = s.find_first_not_of(" \t\r\n");
    if (b == std::string::npos) return {};
    const auto e = s.find_last_not_of(" \t\r\n");
    return s.substr(b, e - b + 1);
}

static std::map<std::string, std::string> read_rsf_header(const fs::path& header_path) {
    std::ifstream in(header_path);
    if (!in) throw std::runtime_error("Cannot open RSF header: " + header_path.string());

    std::ostringstream oss;
    oss << in.rdbuf();
    const std::string text = oss.str();

    // Parse key=value, where value can be quoted or unquoted.
    // This is sufficient for standard Madagascar-style RSF headers.
    std::map<std::string, std::string> kv;
    const std::regex re(R"(([A-Za-z_][A-Za-z0-9_]*)\s*=\s*("[^"]*"|'[^']*'|[^\s]+))");
    for (std::sregex_iterator it(text.begin(), text.end(), re), end; it != end; ++it) {
        std::string key = (*it)[1].str();
        std::string val = trim((*it)[2].str());
        if (val.size() >= 2 && ((val.front() == '"' && val.back() == '"') ||
                                (val.front() == '\'' && val.back() == '\''))) {
            val = val.substr(1, val.size() - 2);
        }
        kv[key] = val;
    }
    return kv;
}

static int get_int(const std::map<std::string, std::string>& h,
                   const std::string& key, int default_value = -1) {
    auto it = h.find(key);
    if (it == h.end()) {
        if (default_value >= 0) return default_value;
        throw std::runtime_error("Missing header key: " + key);
    }
    size_t pos = 0;
    long long v = std::stoll(it->second, &pos);
    if (pos != it->second.size() || v <= 0 || v > std::numeric_limits<int>::max())
        throw std::runtime_error("Invalid integer header value: " + key + "=" + it->second);
    return static_cast<int>(v);
}

static double get_double(const std::map<std::string, std::string>& h,
                         const std::string& key, double default_value) {
    auto it = h.find(key);
    if (it == h.end()) return default_value;
    size_t pos = 0;
    double v = std::stod(it->second, &pos);
    if (pos != it->second.size() || !std::isfinite(v))
        throw std::runtime_error("Invalid floating-point header value: " + key + "=" + it->second);
    return v;
}

static std::string get_string(const std::map<std::string, std::string>& h,
                              const std::string& key,
                              const std::string& default_value = "") {
    auto it = h.find(key);
    return it == h.end() ? default_value : it->second;
}

static std::map<std::string, std::string> parse_args(int argc, char** argv) {
    std::map<std::string, std::string> a;
    for (int i = 1; i < argc; ++i) {
        std::string s(argv[i]);
        const auto p = s.find('=');
        if (p == std::string::npos || p == 0 || p + 1 >= s.size())
            throw std::runtime_error("Arguments must be key=value: " + s);
        a[s.substr(0, p)] = s.substr(p + 1);
    }
    return a;
}

static fs::path resolve_input_binary(const fs::path& header_path,
                                     const std::string& in_value) {
    fs::path p(in_value);

    // 1) Absolute path.
    if (p.is_absolute() && fs::exists(p)) return p;

    // 2) Resolve exactly as written against the current working directory.
    if (fs::exists(p)) return fs::absolute(p);

    // 3) Resolve relative to the RSF header directory.
    fs::path q = header_path.parent_path() / p;
    if (fs::exists(q)) return fs::absolute(q);

    std::ostringstream msg;
    msg << "Cannot find RSF binary file referenced by in=\"" << in_value << "\".\n"
        << "Tried:\n  " << fs::absolute(p).string() << "\n  " << fs::absolute(q).string();
    throw std::runtime_error(msg.str());
}

static fs::path output_binary_path(const fs::path& output_header) {
    return fs::path(output_header.string() + "@");
}

static void ensure_parent(const fs::path& p) {
    if (!p.parent_path().empty()) fs::create_directories(p.parent_path());
}

static std::string quote(const std::string& s) {
    return "\"" + s + "\"";
}

int main(int argc, char** argv) {
    try {
        const auto args = parse_args(argc, argv);
        if (!args.count("input") || !args.count("output") ||
            (!args.count("ix") && !args.count("x"))) {
            std::cerr
                << "Usage:\n"
                << "  " << argv[0] << " input=full_adcig.rsf output=gather.rsf ix=384\n"
                << "or\n"
                << "  " << argv[0] << " input=full_adcig.rsf output=gather.rsf x=3.84\n\n"
                << "Input axes : n1=depth, n2=lateral position, n3=angle\n"
                << "Output axes: n1=depth, n2=angle\n";
            return 1;
        }

        const fs::path input_header = fs::absolute(args.at("input"));
        const fs::path output_header = fs::absolute(args.at("output"));
        const auto h = read_rsf_header(input_header);

        const int n1 = get_int(h, "n1");
        const int n2 = get_int(h, "n2");
        const int n3 = get_int(h, "n3");
        const int esize = get_int(h, "esize", 4);

        const double o1 = get_double(h, "o1", 0.0);
        const double d1 = get_double(h, "d1", 1.0);
        const double o2 = get_double(h, "o2", 0.0);
        const double d2 = get_double(h, "d2", 1.0);
        const double o3 = get_double(h, "o3", 0.0);
        const double d3 = get_double(h, "d3", 1.0);

        const std::string format = get_string(h, "data_format", "native_float");
        if (esize != 4)
            throw std::runtime_error("Only esize=4 float RSF data are supported; got esize=" + std::to_string(esize));
        if (format.find("float") == std::string::npos)
            throw std::runtime_error("Only float RSF data are supported; got data_format=" + format);

        int ix = -1;
        if (args.count("ix")) {
            size_t pos = 0;
            long long tmp = std::stoll(args.at("ix"), &pos);
            if (pos != args.at("ix").size() || tmp < 0 || tmp >= n2)
                throw std::runtime_error("ix must be in [0," + std::to_string(n2 - 1) + "]");
            ix = static_cast<int>(tmp);
        } else {
            if (d2 == 0.0) throw std::runtime_error("Cannot use x= because d2=0");
            size_t pos = 0;
            const double x_req = std::stod(args.at("x"), &pos);
            if (pos != args.at("x").size() || !std::isfinite(x_req))
                throw std::runtime_error("Invalid x= value");
            const long long nearest = std::llround((x_req - o2) / d2);
            if (nearest < 0 || nearest >= n2) {
                std::ostringstream msg;
                msg << "x=" << x_req << " is outside lateral axis ["
                    << o2 << ", " << (o2 + (n2 - 1) * d2) << "]";
                throw std::runtime_error(msg.str());
            }
            ix = static_cast<int>(nearest);
        }
        const double x_actual = o2 + static_cast<double>(ix) * d2;

        const std::string in_value = get_string(h, "in");
        if (in_value.empty()) throw std::runtime_error("Input RSF header has no in= field");
        const fs::path input_binary = resolve_input_binary(input_header, in_value);

        const uint64_t expected_bytes = static_cast<uint64_t>(n1) * n2 * n3 * sizeof(float);
        const uint64_t actual_bytes = fs::file_size(input_binary);
        if (actual_bytes < expected_bytes) {
            std::ostringstream msg;
            msg << "Input binary is too small: " << actual_bytes
                << " bytes, expected at least " << expected_bytes << " bytes";
            throw std::runtime_error(msg.str());
        }

        const fs::path out_binary = output_binary_path(output_header);
        ensure_parent(output_header);
        ensure_parent(out_binary);

        std::ifstream bin_in(input_binary, std::ios::binary);
        if (!bin_in) throw std::runtime_error("Cannot open input binary: " + input_binary.string());
        std::ofstream bin_out(out_binary, std::ios::binary | std::ios::trunc);
        if (!bin_out) throw std::runtime_error("Cannot create output binary: " + out_binary.string());

        // RSF ordering: axis 1 is fastest. For fixed ix and each angle i3,
        // the whole depth trace [0..n1-1] is contiguous.
        std::vector<float> trace(static_cast<size_t>(n1));
        for (int i3 = 0; i3 < n3; ++i3) {
            const uint64_t sample_offset =
                (static_cast<uint64_t>(i3) * static_cast<uint64_t>(n2) +
                 static_cast<uint64_t>(ix)) * static_cast<uint64_t>(n1);
            const uint64_t byte_offset = sample_offset * sizeof(float);

            bin_in.seekg(static_cast<std::streamoff>(byte_offset), std::ios::beg);
            if (!bin_in) throw std::runtime_error("seekg failed while reading angle index " + std::to_string(i3));

            bin_in.read(reinterpret_cast<char*>(trace.data()),
                        static_cast<std::streamsize>(n1 * sizeof(float)));
            if (bin_in.gcount() != static_cast<std::streamsize>(n1 * sizeof(float)))
                throw std::runtime_error("Short read at angle index " + std::to_string(i3));

            bin_out.write(reinterpret_cast<const char*>(trace.data()),
                          static_cast<std::streamsize>(n1 * sizeof(float)));
            if (!bin_out) throw std::runtime_error("Write failed at angle index " + std::to_string(i3));
        }
        bin_out.close();

        // Write a clean 2-D depth-angle RSF header.
        std::ofstream hout(output_header, std::ios::trunc);
        if (!hout) throw std::runtime_error("Cannot create output header: " + output_header.string());

        // Use an absolute binary path for maximum robustness.
        hout << "in=" << quote(out_binary.string()) << "\n";
        hout << "data_format=\"native_float\"\n";
        hout << "esize=4\n";
        hout << std::setprecision(12);
        hout << "n1=" << n1 << "\n";
        hout << "d1=" << d1 << "\n";
        hout << "o1=" << o1 << "\n";
        hout << "label1=" << quote(get_string(h, "label1", "Depth")) << "\n";
        if (!get_string(h, "unit1").empty())
            hout << "unit1=" << quote(get_string(h, "unit1")) << "\n";

        hout << "n2=" << n3 << "\n";
        hout << "d2=" << d3 << "\n";
        hout << "o2=" << o3 << "\n";
        hout << "label2=" << quote(get_string(h, "label3", "Half-opening angle")) << "\n";
        if (!get_string(h, "unit3").empty())
            hout << "unit2=" << quote(get_string(h, "unit3")) << "\n";

        hout << "source_ix=" << ix << "\n";
        hout << "source_x=" << x_actual << "\n";
        hout.close();

        const uint64_t out_expected = static_cast<uint64_t>(n1) * n3 * sizeof(float);
        const uint64_t out_actual = fs::file_size(out_binary);
        if (out_actual != out_expected) {
            std::ostringstream msg;
            msg << "Output binary size mismatch: got " << out_actual
                << ", expected " << out_expected;
            throw std::runtime_error(msg.str());
        }

        std::cout << "ADCIG extraction finished.\n"
                  << "Input volume : " << input_header.string() << "\n"
                  << "Input dims   : n1=" << n1 << " n2=" << n2 << " n3=" << n3 << "\n"
                  << "Selected     : ix=" << ix << ", x=" << x_actual << "\n"
                  << "Output gather: " << output_header.string() << "\n"
                  << "Output dims  : n1=" << n1 << " n2=" << n3 << " (depth x angle)\n"
                  << "Binary bytes : " << out_actual << "\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << "\n";
        return 2;
    }
}
