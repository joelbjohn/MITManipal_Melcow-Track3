#!/usr/bin/env python3
"""IMU-only Track 3 estimator. Usage: python src/main.py INPUT [OUTPUT].

Accepts a EuRoC imu0/data.csv, the organiser's degraded seven-column CSV,
a directory containing imu0/data.csv, or a ZIP containing imu0/data.csv.
Uses no ground truth, waypoint information, camera, or third-party packages.

v2: orientation = gravity-corrected gyro integration (unchanged). Position = offline, zero-phase
"double leaky integrator" of world-frame acceleration (band-limited double integration, cutoff
--fc Hz), because IMU-only integration drifts without bound and only the band above ~0.1 Hz is
recoverable. Falls back to a larger cutoff, then to a constant position, if the result blows up.
"""
import argparse
import csv
import io
import math
import statistics
import zipfile
from array import array
from contextlib import ExitStack
from pathlib import Path

GRAVITY = 9.80665


def multiply(q, r):
    w, x, y, z = q
    a, b, c, d = r
    return (w*a-x*b-y*c-z*d, w*b+x*a+y*d-z*c,
            w*c-x*d+y*a+z*b, w*d+x*c-y*b+z*a)


def rotate(q, v):
    w, x, y, z = q
    vx, vy, vz = v
    tx, ty, tz = 2*(y*vz-z*vy), 2*(z*vx-x*vz), 2*(x*vy-y*vx)
    return (vx+w*tx+y*tz-z*ty,
            vy+w*ty+z*tx-x*tz,
            vz+w*tz+x*ty-y*tx)


def norm(v):
    return math.sqrt(sum(x*x for x in v))


def integrate_rotation(q, omega, dt):
    angle = norm(omega) * dt
    if angle > 1e-10:
        scale = math.sin(0.5*angle) / norm(omega)
        dq = (math.cos(0.5*angle), *(scale*x for x in omega))
    else:
        dq = (1.0, *(0.5*dt*x for x in omega))
    q = multiply(q, dq)
    length = norm(q)
    return tuple(x/length for x in q)


def despike(a, b, c, threshold):
    # Suppress an isolated axis spike only if its two neighbours agree.
    return tuple((a[i]+c[i])*0.5 if abs(b[i]-a[i]) > threshold
                 and abs(b[i]-c[i]) > threshold
                 and abs(a[i]-c[i]) < threshold else b[i] for i in range(3))


def rows_from_file(stream):
    """Yield (header, (timestamp, six sensor values)) for every data row.
    A row whose sensor fields are missing / non-numeric / non-finite is kept with zeros
    (treated downstream as an invalid sample) so the output still has one pose per row."""
    header = None
    for line_number, row in enumerate(csv.reader(stream), 1):
        if not row or not any(x.strip() for x in row):
            continue
        if row[0].lstrip().startswith('#'):
            if len(row) >= 7 and 'timestamp' in row[0].lower():
                header = row[:7]
            continue
        try:
            timestamp = int(row[0]) if row[0].strip().lstrip('+-').isdigit() else float(row[0])
            if not math.isfinite(float(timestamp)):
                continue
        except ValueError:
            if 'timestamp' in row[0].lower():
                header = row[:7]
            continue
        try:
            sensors = tuple(float(value) for value in row[1:7])
            if len(sensors) != 6 or not all(math.isfinite(v) for v in sensors):
                sensors = (0.0,) * 6
        except ValueError:
            sensors = (0.0,) * 6
        yield header, (timestamp, sensors)


def format_from_samples(header, samples):
    first = second = None
    if header:
        name1 = ' '.join(header[1:4]).lower()
        name2 = ' '.join(header[4:7]).lower()
        def kind(s):
            if any(tag in s for tag in ('a_rs', 'accel', 'ax', 'ay', 'az')):
                return 'acc'
            if any(tag in s for tag in ('w_rs', 'gyro', 'gx', 'gy', 'gz')):
                return 'gyro'
            return None
        first, second = kind(name1), kind(name2)
    if first == 'acc' and second == 'gyro':
        accel_first = True
    elif first == 'gyro' and second == 'acc':
        accel_first = False
    else:
        magnitudes = [statistics.median(norm(row[1][offset:offset+3])
                                        for row in samples[:25]) for offset in (0, 3)]
        errors = [min(abs(m-GRAVITY), GRAVITY*abs(m-1.0)) for m in magnitudes]
        if abs(errors[0]-errors[1]) < 0.5:
            raise ValueError('cannot infer sensor column order; add an IMU header')
        accel_first = errors[0] < errors[1]
    offset = 0 if accel_first else 3
    accel_median = statistics.median(norm(row[1][offset:offset+3]) for row in samples[:25])
    accel_scale = GRAVITY if 0.5 < accel_median < 1.5 else 1.0
    diffs = [b[0]-a[0] for a, b in zip(samples, samples[1:]) if b[0] > a[0]]
    if not diffs:
        raise ValueError('need at least two distinct timestamps')
    raw_dt = statistics.median(diffs)
    ts_scale = 1e-9 if raw_dt > 1e5 else 1e-6 if raw_dt > 1e2 else 1e-3 if raw_dt > 0.5 else 1.0
    if header and '[ns]' in header[0].lower():
        ts_scale = 1e-9
    return accel_first, accel_scale, ts_scale


def timestamp_text(raw, scale):
    # Keep integer nanosecond timestamps exact in decimal seconds.
    if scale == 1e-9 and isinstance(raw, int):
        seconds, nanos = divmod(raw, 1_000_000_000)
        return f'{seconds}.{nanos:09d}'
    return f'{raw*scale:.9f}'


def leaky_pass(buf, ts, ts_scale, nominal_dt, wc, forward):
    """In-place first-order leaky integrator on interleaved xyz data (variable dt).
    Run forward then backward it is zero-phase. Rows with dt == 0 (duplicates) change nothing."""
    n = len(ts)
    exp = math.exp
    y0 = y1 = y2 = 0.0
    order = range(n) if forward else range(n - 1, -1, -1)
    for i in order:
        if forward:
            dt = (ts[i] - ts[i - 1]) * ts_scale if i else nominal_dt
        else:
            dt = (ts[i + 1] - ts[i]) * ts_scale if i < n - 1 else nominal_dt
        e = exp(-wc * dt)
        j = 3 * i
        y0 = e * y0 + dt * buf[j]
        y1 = e * y1 + dt * buf[j + 1]
        y2 = e * y2 + dt * buf[j + 2]
        buf[j] = y0
        buf[j + 1] = y1
        buf[j + 2] = y2


def band_positions(acc_w, ts, ts_scale, nominal_dt, fc, gain, spread_max=4.0):
    """Double leaky integration (4 passes, zero phase). Tries fc, then larger cutoffs if the
    estimate is implausibly large; the last resort is a constant position."""
    n = len(ts)
    for factor in (1.0, 1.5, 2.5, 4.0):
        wc = 2.0 * math.pi * fc * factor
        buf = array('f', acc_w)
        for forward in (True, False, True, False):
            leaky_pass(buf, ts, ts_scale, nominal_dt, wc, forward)
        mean = [sum(buf[k::3]) / n for k in range(3)]
        sample = range(0, n, 4)
        spread = math.sqrt(sum(sum((buf[3*i+k] - mean[k]) ** 2 for k in range(3))
                               for i in sample) / len(sample))
        if spread * gain <= spread_max:
            if gain != 1.0:
                for i in range(len(buf)):
                    buf[i] *= gain
            return buf
    return array('f', bytes(4 * 3 * n))


def solve(source, destination, fc=0.15, gain=1.0):
    if fc <= 0 or gain < 0:
        raise ValueError('fc must be positive and gain nonnegative')
    with ExitStack() as stack:
        path = Path(source)
        if path.is_dir():
            path = path / 'imu0' / 'data.csv'
        if zipfile.is_zipfile(path):
            archive = stack.enter_context(zipfile.ZipFile(path))
            matches = [name for name in archive.namelist()
                       if name.replace('\\', '/').endswith('imu0/data.csv')]
            if len(matches) != 1:
                raise ValueError('ZIP must contain exactly one imu0/data.csv')
            stream = stack.enter_context(io.TextIOWrapper(
                stack.enter_context(archive.open(matches[0])), encoding='utf-8-sig', newline=''))
        else:
            stream = stack.enter_context(open(path, encoding='utf-8-sig', newline=''))
        rows = rows_from_file(stream)
        warmup = []
        header = None
        for h, row in rows:
            if h:
                header = h
            warmup.append(row)
            if len(warmup) >= 20:
                _, _, scale = format_from_samples(header, warmup[:25])
                if (warmup[-1][0]-warmup[0][0])*scale >= 2.0 or len(warmup) >= 450:
                    break
        if len(warmup) < 2:
            raise ValueError('input must contain at least two IMU rows')
        accel_first, accel_scale, ts_scale = format_from_samples(header, warmup[:25])
        times = [b[0]-a[0] for a, b in zip(warmup, warmup[1:]) if b[0] > a[0]]
        nominal_dt = statistics.median(times)*ts_scale

        def parse(row):
            raw, six = row
            if accel_first:
                acc = tuple(v*accel_scale for v in six[:3])
                gyro = six[3:]
            else:
                gyro = six[:3]
                acc = tuple(v*accel_scale for v in six[3:])
            return raw, acc, gyro

        # Average a short initial interval. Only use the full two-second window
        # if both vector acceleration and angular velocity are stable.
        initial = [parse(r) for r in warmup]
        gravity_candidates = [r for r in initial if (r[0]-initial[0][0])*ts_scale <= 0.15]
        def stationary_window(window):
            mean_w = tuple(statistics.fmean(r[2][j] for r in window) for j in range(3))
            mean_a = tuple(statistics.fmean(r[1][j] for r in window) for j in range(3))
            stable = (len(window) >= 10 and abs(norm(mean_a)-GRAVITY) < 0.3
                      and max(norm(tuple(r[2][j]-mean_w[j] for j in range(3)))
                              for r in window) < 0.03
                      and max(norm(tuple(r[1][j]-mean_a[j] for j in range(3)))
                              for r in window) < 0.18)
            return stable, mean_w
        stable, mean_gyr = stationary_window(initial)
        if stable:
            gravity_candidates = initial
        else:
            stable, mean_gyr = stationary_window(gravity_candidates)
        ax, ay, az = (statistics.fmean(r[1][j] for r in gravity_candidates) for j in range(3))
        roll = math.atan2(ay, az)
        pitch = math.atan2(-ax, math.hypot(ay, az))
        cr, sr = math.cos(roll/2), math.sin(roll/2)
        cp, sp = math.cos(pitch/2), math.sin(pitch/2)
        q = (cr*cp, sr*cp, cr*sp, -sr*sp)  # world yaw is arbitrary
        bias = list(mean_gyr) if stable else [0.0, 0.0, 0.0]
        previous_omega = (0.0, 0.0, 0.0)
        kp, ki = 0.4, 0.02

        def all_rows():
            for item in warmup:
                yield parse(item)
            for _, item in rows:
                yield parse(item)

        ts_arr = array('q' if isinstance(initial[0][0], int) else 'd')
        q_arr = array('f')
        a_arr = array('f')

        def record(raw, aw):
            ts_arr.append(raw)
            q_arr.extend(q)
            a_arr.extend(aw)

        iterator = iter(all_rows())
        prev = next(iterator)
        record(prev[0], (0.0, 0.0, 0.0))
        curr = next(iterator, None)
        while curr is not None:
            nxt = next(iterator, None)
            raw, acc_raw, gyro_raw = curr
            dt = (raw-prev[0])*ts_scale
            if dt < 0:
                raise ValueError('input timestamps must be nondecreasing')
            if dt == 0:
                record(raw, (0.0, 0.0, 0.0))
                prev, curr = curr, nxt
                continue

            dead = acc_raw == (0.0, 0.0, 0.0) and gyro_raw == (0.0, 0.0, 0.0)
            stuck = (acc_raw == prev[1] and gyro_raw == prev[2])
            gap = dt > 1.5*nominal_dt
            aw = (0.0, 0.0, 0.0)
            if dead or gap:
                if not dead:
                    # Missing measurements: extrapolate the last rate briefly (capped).
                    q = integrate_rotation(q, previous_omega, min(dt, 0.25))
                    if dt <= 0.3:   # short gap: hold the sample at the end of the gap
                        w = rotate(q, acc_raw)
                        aw = (w[0], w[1], w[2]-GRAVITY)
            else:
                if (nxt is not None and 0 < (nxt[0]-raw)*ts_scale <= 1.5*nominal_dt
                        and prev[1] != (0.0, 0.0, 0.0)
                        and nxt[1] != (0.0, 0.0, 0.0)):
                    acc = despike(prev[1], acc_raw, nxt[1], 3.0)
                    gyro = despike(prev[2], gyro_raw, nxt[2], 0.5)
                else:
                    acc, gyro = acc_raw, gyro_raw
                omega = tuple(gyro[j]-bias[j] for j in range(3))
                magnitude = norm(acc)
                expected = rotate((q[0], -q[1], -q[2], -q[3]), (0., 0., 1.))
                gravity_residual = norm(tuple(acc[j]-GRAVITY*expected[j] for j in range(3)))
                # Norm alone cannot reject horizontal acceleration: it can
                # remain close to g while rotating the gravity direction.
                if (not stuck and abs(magnitude-GRAVITY) < 0.08*GRAVITY
                        and gravity_residual < 0.12 and norm(omega) < 1.0):
                    measured = tuple(x/magnitude for x in acc)
                    error = (measured[1]*expected[2]-measured[2]*expected[1],
                             measured[2]*expected[0]-measured[0]*expected[2],
                             measured[0]*expected[1]-measured[1]*expected[0])
                    for j in range(3):
                        bias[j] -= ki*error[j]*dt
                    omega = tuple(omega[j]+kp*error[j] for j in range(3))
                q = integrate_rotation(q, omega, dt)
                previous_omega = omega
                w = rotate(q, acc)
                aw = (w[0], w[1], w[2]-GRAVITY)
            record(raw, aw)
            prev, curr = curr, nxt

        pos = band_positions(a_arr, ts_arr, ts_scale, nominal_dt, fc, gain)
        with open(destination, 'w', encoding='utf-8') as out:
            out.write('# timestamp tx ty tz qx qy qz qw\n')
            for i in range(len(ts_arr)):
                j, k = 3*i, 4*i
                out.write(f'{timestamp_text(ts_arr[i], ts_scale)} {pos[j]:.8f} {pos[j+1]:.8f} {pos[j+2]:.8f} '
                          f'{q_arr[k+1]:.9f} {q_arr[k+2]:.9f} {q_arr[k+3]:.9f} {q_arr[k]:.9f}\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', help='CSV, ZIP, or directory with imu0/data.csv')
    parser.add_argument('output', nargs='?', default='estimated_trajectory.txt')
    parser.add_argument('--fc', type=float, default=0.15,
                        help='cutoff of the double leaky integrator in Hz; default 0.15')
    parser.add_argument('--gain', type=float, default=1.0, help='position gain; default 1')
    args = parser.parse_args()
    solve(args.input, args.output, args.fc, args.gain)


if __name__ == '__main__':
    main()
