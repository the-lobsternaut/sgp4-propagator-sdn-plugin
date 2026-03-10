/**
 * SGP4 Propagator Implementation
 *
 * Uses dnwrnr/sgp4 (Apache 2.0) via direct OrbitalElements constructor.
 * Bypasses TLE text entirely — integer NORAD IDs with no format limit.
 */

#include "sgp4_prop/propagator.h"

// dnwrnr/sgp4 headers
#include "OrbitalElements.h"
#include "SGP4.h"
#include "DateTime.h"
#include "TimeSpan.h"
#include "Eci.h"

#include <cmath>
#include <stdexcept>

namespace sgp4_prop {

static constexpr double DEG_TO_RAD = 0.017453292519943295;
static constexpr double TWOPI = 6.283185307179586;
static constexpr double SEC_PER_DAY = 86400.0;
static constexpr double MIN_PER_DAY = 1440.0;
static constexpr double XKE = 0.0743669161;  // sqrt(GM) in SGP4 units

// ── Julian Date ↔ DateTime conversion ──

static libsgp4::DateTime jd_to_datetime(double jd) {
    // JD → calendar date
    double jd_plus = jd + 0.5;
    int z = static_cast<int>(jd_plus);
    double f = jd_plus - z;

    int a;
    if (z < 2299161) { a = z; }
    else {
        int alpha = static_cast<int>((z - 1867216.25) / 36524.25);
        a = z + 1 + alpha - alpha / 4;
    }

    int b = a + 1524;
    int c = static_cast<int>((b - 122.1) / 365.25);
    int d = static_cast<int>(365.25 * c);
    int e = static_cast<int>((b - d) / 30.6001);

    int day = b - d - static_cast<int>(30.6001 * e);
    int month = (e < 14) ? e - 1 : e - 13;
    int year = (month > 2) ? c - 4716 : c - 4715;

    double day_frac = f;
    int hours = static_cast<int>(day_frac * 24.0);
    double rem = day_frac * 24.0 - hours;
    int minutes = static_cast<int>(rem * 60.0);
    rem = rem * 60.0 - minutes;
    int seconds = static_cast<int>(rem * 60.0);
    int microseconds = static_cast<int>((rem * 60.0 - seconds) * 1000000);

    return libsgp4::DateTime(year, month, day, hours, minutes, seconds, microseconds);
}

// ── Direct GP → SGP4 propagation (no TLE text) ──

StateVector propagate_to_epoch(const GPElement& gp, double target_jd) {
    // Convert GP to OrbitalElements directly (8-arg constructor)
    double n_rad_min = gp.mean_motion * TWOPI / MIN_PER_DAY;

    libsgp4::DateTime epoch_dt = jd_to_datetime(gp.epoch_jd);

    libsgp4::OrbitalElements elements(
        gp.mean_anomaly * DEG_TO_RAD,
        gp.ra_of_asc_node * DEG_TO_RAD,
        gp.arg_of_pericenter * DEG_TO_RAD,
        gp.eccentricity,
        gp.inclination * DEG_TO_RAD,
        n_rad_min,
        gp.bstar,
        epoch_dt
    );

    libsgp4::SGP4 sgp4(elements);

    // Use DateTime for FindPosition
    libsgp4::DateTime target_dt = jd_to_datetime(target_jd);

    libsgp4::Eci eci = sgp4.FindPosition(target_dt);
    auto pos = eci.Position();
    auto vel = eci.Velocity();

    StateVector sv;
    sv.epoch_jd = target_jd;
    sv.x = pos.x; sv.y = pos.y; sv.z = pos.z;
    sv.vx = vel.x; sv.vy = vel.y; sv.vz = vel.z;
    return sv;
}

// ── Range propagation ──

std::vector<StateVector> propagate(const GPElement& gp, const PropagationConfig& config) {
    std::vector<StateVector> states;

    double start_jd = config.start_jd;
    if (start_jd <= 0.0) start_jd = gp.epoch_jd;

    double end_jd = config.end_jd;
    if (end_jd <= 0.0) end_jd = start_jd + config.duration_days;

    double step_days = config.step_seconds / SEC_PER_DAY;
    int n_steps = static_cast<int>((end_jd - start_jd) / step_days) + 1;
    states.reserve(n_steps);

    // Set up SGP4 once
    double n_rad_min = gp.mean_motion * TWOPI / MIN_PER_DAY;

    libsgp4::DateTime epoch_dt = jd_to_datetime(gp.epoch_jd);
    libsgp4::OrbitalElements elements(
        gp.mean_anomaly * DEG_TO_RAD, gp.ra_of_asc_node * DEG_TO_RAD,
        gp.arg_of_pericenter * DEG_TO_RAD, gp.eccentricity,
        gp.inclination * DEG_TO_RAD, n_rad_min, gp.bstar, epoch_dt
    );
    libsgp4::SGP4 sgp4(elements);

    for (double t = start_jd; t <= end_jd; t += step_days) {
        libsgp4::DateTime target_dt = jd_to_datetime(t);

        try {
            libsgp4::Eci eci = sgp4.FindPosition(target_dt);
            auto pos = eci.Position();
            auto vel = eci.Velocity();

            StateVector sv;
            sv.epoch_jd = t;
            sv.x = pos.x; sv.y = pos.y; sv.z = pos.z;
            sv.vx = vel.x; sv.vy = vel.y; sv.vz = vel.z;
            states.push_back(sv);
        } catch (...) {
            // SGP4 can fail for decayed objects — skip this step
        }
    }

    return states;
}

// ── OEM output (stub — will use FlatBuffers when wired up) ──

int32_t propagate_to_oem(const GPElement& gp,
                          const PropagationConfig& config,
                          uint8_t* output, uint32_t output_capacity) {
    auto states = propagate(gp, config);
    if (states.empty()) return -3;  // invalid input

    // TODO: Build OEM FlatBuffers from states
    // For now, return state count as success indicator
    return static_cast<int32_t>(states.size());
}

int32_t propagate_batch_to_oem(const std::vector<GPElement>& gps,
                                const PropagationConfig& config,
                                uint8_t* output, uint32_t output_capacity) {
    // TODO: Build OEM collection FlatBuffers
    int32_t total = 0;
    for (const auto& gp : gps) {
        auto states = propagate(gp, config);
        total += static_cast<int32_t>(states.size());
    }
    return total;
}

}  // namespace sgp4_prop
