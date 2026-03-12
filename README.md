# SGP4 Propagator SDN Plugin

SGP4/SDP4 orbit propagator for the [Space Data Network](https://github.com/the-lobsternaut/space-data-network). Takes GP/OMM elements (CelesTrak JSON format), propagates with the standard AFSPC analytical theory, and outputs OEM FlatBuffers.

## Features

- **Two propagation paths**:
  - GP JSON → direct `OrbitalElements` → SGP4 (integer NORAD IDs, no 5-digit limit)
  - TLE text → SGP4 (legacy, Alpha-5 up to 339999)
- **Batch propagation** — process entire GP catalogs in one call
- **Point propagation** — single-epoch queries for CA/OD plugins
- **OEM output** — spacedatastandards.org FlatBuffers (`$OEM` file identifier)
- **WASM-ready** — Emscripten compilation for browser/Node.js

## Usage

```cpp
#include "sgp4_prop/propagator.h"

// Parse CelesTrak GP JSON
auto elements = sgp4_prop::parse_gp_json(json_string);

// Propagate to OEM
sgp4_prop::PropagationConfig config;
config.duration_days = 7.0;
config.step_seconds = 60.0;

auto states = sgp4_prop::propagate(elements[0], config);
// states: vector of StateVector (epoch_jd, x, y, z, vx, vy, vz in TEME km/s)

// Point propagation (fast single-epoch query)
auto state = sgp4_prop::propagate_to_epoch(elements[0], target_jd);

// OEM FlatBuffer output
uint8_t buffer[1024*1024];
int32_t size = sgp4_prop::propagate_to_oem(elements[0], config, buffer, sizeof(buffer));
// buffer contains $OEM FlatBuffer binary
```

## Data Flow

```
CelesTrak GP/OMM JSON ──→ parse_gp_json() ──→ GPElement
                                                  │
                                    propagate() / propagate_to_epoch()
                                                  │
                                                  ▼
                                          StateVector (TEME)
                                                  │
                                      propagate_to_oem()
                                                  │
                                                  ▼
                                        $OEM FlatBuffer binary
```

## GP Element Fields

| Field | Type | Description |
|---|---|---|
| `norad_cat_id` | int | NORAD catalog ID (integer, no format limit) |
| `mean_motion` | double | Mean motion (rev/day) |
| `eccentricity` | double | Eccentricity |
| `inclination` | double | Inclination (degrees) |
| `ra_of_asc_node` | double | RAAN (degrees) |
| `arg_of_pericenter` | double | Argument of perigee (degrees) |
| `mean_anomaly` | double | Mean anomaly (degrees) |
| `bstar` | double | B* drag term |
| `epoch` | string | ISO 8601 epoch |

## Building

```bash
cd src/cpp && mkdir -p build && cd build
cmake ..
make -j4
./test_propagator
```

## Dependencies

- [dnwrnr/sgp4](https://github.com/dnwrnr/sgp4) — C++ SGP4/SDP4 implementation (Apache-2.0)
- [spacedatastandards.org](https://spacedatastandards.org) — OEM FlatBuffers schema

## Coordinate System

Output states are in **TEME** (True Equator Mean Equinox), the native frame for SGP4. Units: km and km/s.

## License

Apache-2.0
