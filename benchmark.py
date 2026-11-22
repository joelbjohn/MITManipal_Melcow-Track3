#!/usr/bin/env python3
"""Local Track 3 budget check. Run from repository root:
python3 benchmark.py data/degraded_data.csv
Writes estimated_trajectory.txt and metadata_local.json.
"""
import argparse
import csv
import importlib.util
import json
import math
import statistics
import time
import tracemalloc
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='degraded seven-column IMU CSV')
    parser.add_argument('--output', default='estimated_trajectory.txt')
    parser.add_argument('--metadata', default='metadata.json')
    parser.add_argument('--runs', type=int, default=3)
    args = parser.parse_args()
    if args.runs < 1:
        parser.error('--runs must be positive')
    script = Path(__file__).resolve().parent / 'src' / 'main.py'
    if not script.is_file():
        parser.error(f'missing {script}; save benchmark.py in the repository root')
    spec = importlib.util.spec_from_file_location('track3_main', script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    with open(args.input, newline='', encoding='utf-8-sig') as file:
        rows = csv.reader(file)
        times = []
        for row in rows:
            if not row or row[0].startswith('#') or 'timestamp' in row[0].lower():
                continue
            times.append(row[0].strip())
    if len(times) < 2:
        raise ValueError('need two or more IMU rows')
    raw_dt = float(times[1]) - float(times[0])
    scale = 1e-9 if raw_dt > 1e5 else 1e-6 if raw_dt > 1e2 else 1e-3 if raw_dt > 0.5 else 1.0
    # Integer subtraction preserves precision for epoch nanoseconds.
    diff = int(times[-1]) - int(times[0]) if scale == 1e-9 else float(times[-1]) - float(times[0])
    duration = diff * scale
    if duration <= 0:
        raise ValueError('IMU duration must be positive')

    runtimes = []
    for _ in range(args.runs):
        begin = time.perf_counter()
        module.solve(args.input, args.output)
        runtimes.append(time.perf_counter() - begin)
    tracemalloc.start()
    module.solve(args.input, args.output)
    peak_bytes = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()

    output_count = 0
    with open(args.output, encoding='utf-8') as file:
        for line in file:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) != 8 or not all(math.isfinite(float(v)) for v in fields):
                raise ValueError(f'invalid output pose at row {output_count + 1}')
            output_count += 1
    runtime = statistics.median(runtimes)
    budget_seconds = duration * 0.020
    result = {
        'input': str(Path(args.input).resolve()),
        'imu_duration_s': round(duration, 6),
        'rows_input': len(times),
        'rows_output': output_count,
        'output_rows_match': output_count == len(times),
        'runtime_s_median': round(runtime, 6),
        'runtime_ms_per_imu_s': round(runtime * 1000 / duration, 4),
        'runtime_budget_s': round(budget_seconds, 6),
        'runtime_budget_pass_local': runtime <= budget_seconds,
        'peak_traced_python_allocations_bytes': peak_bytes,
        'heap_budget_bytes': 8 * 1024 * 1024,
        'heap_budget_pass_local': peak_bytes <= 8 * 1024 * 1024,
        'note': 'Local checks only. tracemalloc excludes some native allocations and base interpreter RSS; organiser scoring may differ.',
    }
    with open(args.metadata, 'w', encoding='utf-8') as file:
        json.dump(result, file, indent=2)
    print(json.dumps(result, indent=2))
    if not result['output_rows_match']:
        raise SystemExit('FAIL: output row count differs from input')


if __name__ == '__main__':
    main()
