#!/usr/bin/env bash
# Run the calibration-logic tests (no susieR needed) then the full benchmark.
set -e
echo "=== Calibration logic tests ==="
Rscript tests/test_calibration.R
echo ""
echo "=== Full benchmark (requires susieR) ==="
Rscript R/run_benchmark.R
