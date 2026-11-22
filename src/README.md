# Path Finder: Estimating a Flight Path from a Damaged IMU

**Samsung XR Hackathon 2026, Track 3**  
**Team:** MITManipal_Melcow-Track3

## What we built

We wrote a Python program that estimates where a flying device moved using
only its IMU readings: acceleration and rotation over time. The challenge
damages these readings as time passes. Some readings become noisy or repeat,
some rows disappear, and the sensor values turn to zero for the last five
seconds. Our program still writes a position and orientation for every row
it receives.

The program does not read the known path or any checkpoint locations. We use
the known path from the public EuRoC dataset only afterwards to see how far
our estimate is from it.

## Run it

From the project folder, run:

```bash
python3 src/main.py data/degraded_data.csv estimated_trajectory.txt
```

You can replace the input and output names with other paths. `main.py` uses
Python 3's standard library, so it does not need a model download or extra
Python packages. It accepts a seven-column IMU CSV, a ZIP containing
`imu0/data.csv`, or a folder containing `imu0/data.csv`.

The output `estimated_trajectory.txt` has one space-separated row for each
input row:

```text
timestamp_seconds x y z qx qy qz qw
```

The input EuRoC timestamps are in nanoseconds. `main.py` writes them in
seconds. A quaternion (`qx qy qz qw`) describes orientation. The starting
position and compass heading are arbitrary because the IMU alone does not
provide a global map. We account for this when comparing paths.

## How it works

The program first identifies the IMU columns and estimates the initial tilt.
If the first readings are steady enough, it also estimates the gyro's initial
offset. It removes isolated spikes when the readings before and after a spike
agree. It then updates orientation from the gyro, using the accelerometer for
small tilt corrections when the reading looks reliable.

Next, it converts acceleration into a fixed world direction and removes
gravity. We process the full recording in a compact buffer and filter out
slow changes that could otherwise make the estimated position drift far
away. This filter has a cutoff of `0.15 Hz` by default. The program predicts
through missing timestamp intervals, repeats the pose for duplicate
timestamps, and estimates the final five seconds when the IMU is zeroed.

This is a tradeoff: filtering drift can also remove real, slow movement.

## Data and local checks

We used two public EuRoC MAV recordings: `V1_01_easy` and `V1_02_medium`.
For each recording, we rearranged a copy of the IMU columns because the
organiser's `degrade_imu.py` expects acceleration before gyro readings.
The organiser's script then produced damaged IMU data with increasing noise,
drift, spikes, repeated values, dropped rows, and a final five seconds of
zero sensor values. We ran `main.py` on that damaged file only. We compared
its output with EuRoC ground truth after the run, aligning the two paths by
rotation and translation, with no scale correction.

| Public recording | Damaged IMU rows | Output poses | Position error (RMSE) | Final 5 s error (RMSE) |
| --- | ---: | ---: | ---: | ---: |
| `V1_01_easy` | 28,367 | 28,367 | 1.75 m | 1.44 m |
| `V1_02_medium` | 16,624 | 16,624 | 1.59 m | 0.73 m |

These are **local public-dataset results**, not scores from the organiser's
hidden checkpoints. The known EuRoC path is used only by the comparison,
never by `main.py`.

On our laptop, the `V1_02_medium` run took a median of **0.249 s** for
**85.495 s of IMU data**, or **2.92 ms of processing per second of IMU**.
The challenge limit is 20 ms per IMU second. Python allocation tracing
reported a peak of **1,127,503 bytes (about 1.08 MiB)**, below the 8 MB
algorithm-memory limit. The Python interpreter itself is excluded by the
organiser's FAQ. These measurements are from `benchmark.py`; the organiser
will measure performance in its own environment.

To measure this again on your computer, use:

```bash
python3 benchmark.py data/degraded_data.csv
```

This creates `metadata_local.json` and regenerates the output trajectory.
`benchmark.py` is for local checking; the evaluator runs `src/main.py`.

## Challenge tiers and limits

- **Checkpoints:** The program produces a continuous path, but the hidden
  checkpoint locations are unavailable to us. We cannot report a checkpoint
  score.
- **Final five seconds:** The program writes a pose for every zeroed IMU row.
  We measured its prediction error on both public recordings above.
- **Dropouts:** Missing rows appear as gaps between timestamps. The program
  detects these gaps and predicts through them. We have not measured a
  separate dropout-only score.

The position errors should be read alongside the path plots. The filter
keeps the path from drifting without bound, but its estimated path misses
some real movement. An IMU by itself cannot tell us the exact starting
speed or completely separate slow movement from sensor error. We have not
claimed a score on the organisers' private data or measured rotation error
separately.

## Repository and submission files

The repository can keep `benchmark.py`, the organiser's degradation script,
and instructions explaining how the sample data was made. Keep the large
EuRoC CSVs and `.venv/` out of GitHub. The final submission ZIP follows the
Track 3 layout:

```text
estimated_trajectory.txt
metadata.json                 #local measurement record
src/main.py
src/README.md               
```

The final video and presentation accompany the submission separately. The
team ZIP and presentation names should match the registered team name.
