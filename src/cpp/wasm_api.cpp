/**
 * sgp4-propagator-sdn-plugin WASM API
 *
 * SGP4/SDP4 propagation from GP/OMM JSON, Lambert STM covariance.
 * JSON-in/JSON-out for propagation; C-API for binary OEM output.
 */

#include "sgp4_prop/propagator.h"
#include "sgp4_prop/lambert_stm.h"

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#include <emscripten/bind.h>
using namespace emscripten;
#endif

#include <string>
#include <cstring>
#include <sstream>
#include <cstdlib>

// ============================================================================
// Version
// ============================================================================

static std::string version() {
    return "0.1.0";
}

// ============================================================================
// JSON-in/JSON-out wrappers
// ============================================================================

/// Parse GP JSON and propagate, returning state vectors as JSON array
static std::string propagate_gp_json(const std::string& gp_json,
                                      double duration_days,
                                      double step_seconds) {
    try {
        auto gps = sgp4_prop::parse_gp_json(gp_json);
        if (gps.empty()) return "{\"error\":\"no GP elements parsed\"}";

        sgp4_prop::PropagationConfig config;
        config.duration_days = duration_days;
        config.step_seconds = step_seconds;

        auto states = sgp4_prop::propagate(gps[0], config);

        std::ostringstream oss;
        oss.precision(15);
        oss << "{\"object_name\":\"" << gps[0].object_name
            << "\",\"norad_id\":" << gps[0].norad_cat_id
            << ",\"num_states\":" << states.size()
            << ",\"states\":[";
        for (size_t i = 0; i < states.size(); ++i) {
            if (i > 0) oss << ",";
            oss << "{\"jd\":" << states[i].epoch_jd
                << ",\"x\":" << states[i].x << ",\"y\":" << states[i].y
                << ",\"z\":" << states[i].z
                << ",\"vx\":" << states[i].vx << ",\"vy\":" << states[i].vy
                << ",\"vz\":" << states[i].vz << "}";
        }
        oss << "]}";
        return oss.str();
    } catch (...) {
        return "{\"error\":\"propagation failed\"}";
    }
}

/// Propagate to single epoch from GP JSON
static std::string propagate_to_epoch_json(const std::string& gp_json,
                                            double target_jd) {
    try {
        auto gps = sgp4_prop::parse_gp_json(gp_json);
        if (gps.empty()) return "{\"error\":\"no GP elements parsed\"}";

        auto sv = sgp4_prop::propagate_to_epoch(gps[0], target_jd);

        std::ostringstream oss;
        oss.precision(15);
        oss << "{\"epoch_jd\":" << sv.epoch_jd
            << ",\"x\":" << sv.x << ",\"y\":" << sv.y << ",\"z\":" << sv.z
            << ",\"vx\":" << sv.vx << ",\"vy\":" << sv.vy << ",\"vz\":" << sv.vz
            << "}";
        return oss.str();
    } catch (...) {
        return "{\"error\":\"propagation failed\"}";
    }
}

/// Solve Lambert problem — direct parameters
static std::string solve_lambert_json(double r1x, double r1y, double r1z,
                                       double r2x, double r2y, double r2z,
                                       double tof_seconds) {
    double r1[3] = {r1x, r1y, r1z};
    double r2[3] = {r2x, r2y, r2z};

    auto sol = sgp4_prop::solve_lambert(r1, r2, tof_seconds);

    std::ostringstream oss;
    oss.precision(15);
    oss << "{\"converged\":" << (sol.converged ? "true" : "false")
        << ",\"iterations\":" << sol.iterations
        << ",\"v1\":[" << sol.v1[0] << "," << sol.v1[1] << "," << sol.v1[2] << "]"
        << ",\"v2\":[" << sol.v2[0] << "," << sol.v2[1] << "," << sol.v2[2] << "]"
        << "}";
    return oss.str();
}

/// Epoch string to JD
static double epoch_to_jd(const std::string& epoch_str) {
    return sgp4_prop::epoch_to_jd(epoch_str);
}

// ============================================================================
// SDN Plugin ABI (C exports)
// ============================================================================

#ifdef __EMSCRIPTEN__

extern "C" {

EMSCRIPTEN_KEEPALIVE
void* sdn_malloc(size_t size) {
    return malloc(size);
}

EMSCRIPTEN_KEEPALIVE
void sdn_free(void* ptr) {
    free(ptr);
}

/// Parse GP JSON + propagate → OEM FlatBuffers binary
EMSCRIPTEN_KEEPALIVE
int32_t propagate_oem(const char* gp_json, size_t json_len,
                      double duration_days, double step_s,
                      uint8_t* output, size_t output_len) {
    try {
        std::string json(gp_json, json_len);
        auto gps = sgp4_prop::parse_gp_json(json);
        if (gps.empty()) return -3;

        sgp4_prop::PropagationConfig config;
        config.duration_days = duration_days;
        config.step_seconds = step_s;

        return sgp4_prop::propagate_to_oem(
            gps[0], config, output, static_cast<uint32_t>(output_len));
    } catch (...) {
        return -1;
    }
}

/// Batch propagate → OEM collection
EMSCRIPTEN_KEEPALIVE
int32_t propagate_batch_oem(const char* gp_json, size_t json_len,
                            double duration_days, double step_s,
                            uint8_t* output, size_t output_len) {
    try {
        std::string json(gp_json, json_len);
        auto gps = sgp4_prop::parse_gp_json(json);
        if (gps.empty()) return -3;

        sgp4_prop::PropagationConfig config;
        config.duration_days = duration_days;
        config.step_seconds = step_s;

        return sgp4_prop::propagate_batch_to_oem(
            gps, config, output, static_cast<uint32_t>(output_len));
    } catch (...) {
        return -1;
    }
}

} // extern "C"

// ============================================================================
// Embind Exports
// ============================================================================

EMSCRIPTEN_BINDINGS(sdn_sgp4_propagator) {
    function("version", &version);
    function("propagateGP", &propagate_gp_json);
    function("propagateToEpoch", &propagate_to_epoch_json);
    function("solveLambert", &solve_lambert_json);
    function("epochToJd", &epoch_to_jd);
}

#endif
