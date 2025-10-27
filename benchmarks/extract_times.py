import os
import sys


def extract_times(file_path, start_str):
    matches = []
    try:
        with open(file_path, "r") as f:
            for line in f:
                if start_str in line:
                    start_index = line.find(start_str) + len(start_str)
                    end_index = line.find("s", start_index)
                    if end_index == -1:
                        continue
                    substring = line[start_index:end_index]
                    try:
                        valueAtTime = float(substring.strip().replace(",", "."))
                        matches.append(valueAtTime)
                    except ValueError:
                        # ignore malformed lines but print a warning
                        print(f"CANNOT CONVERT SUBSTRING to float in {file_path}: '{substring.strip()}'")
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
    return matches

if len(sys.argv) != 2:
    print(
        "Usage: python extract_times.py <folder_path>\n"
        + "<folder_path> the path to the folder containing the log files\n"
    )
    sys.exit(1)

folder_path = sys.argv[1]

files = os.listdir(folder_path)
start_str = "INFO: compute() took"

# Collect first (full) and second (short) occurrence per file
times_full = []  # first occurrence per file
times = []       # second occurrence per file

for file_name in files:
    file_path = os.path.join(folder_path, file_name)
    # skip directories
    if not os.path.isfile(file_path):
        continue
    matches = extract_times(file_path, start_str)
    if len(matches) >= 1:
        times_full.append(matches[0])
    if len(matches) >= 2:
        times.append(matches[1])
    if len(matches) > 2:
        print(f"WARNING: {file_name} contains {len(matches)} matches, only first two used")

if times_full:
    with open("times_full.txt", "w") as file:
        average_full = sum(times_full) / len(times_full)
        file.write(f"Average: {average_full}s\n")
        for i, t in enumerate(times_full):
            file.write(f"Run {i}: {t}s\n")

if times:
    with open("times.txt", "w") as file:
        average = sum(times) / len(times)
        file.write(f"Average: {average}s\n")
        for i, t in enumerate(times):
            file.write(f"Run {i}: {t}s\n")

if len(times_full) != len(files):
    print(f"WARNING: Extracted {len(times_full)} first-times, but there are {len(files)} files in total!")
print(f"Collected {len(times_full)} first-times and {len(times)} second-times from {len(files)} files.")
