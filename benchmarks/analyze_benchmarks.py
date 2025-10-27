import sys
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatter


# do not include header in startRow!
def readInResults(pathToSheet, startRow, endRow, colRange):
    nrows = endRow - startRow + 1
    df = pd.read_excel(
        pathToSheet,
        engine="odf",
        usecols=colRange,
        header=None,
        skiprows=startRow - 1,
        nrows=nrows,
    )

    arrays = []
    for i in range(8):
        col = pd.to_numeric(df.iloc[:, i], errors="coerce")
        if i < 2:
            # integer columns: fail if there are missing/non-finite values
            if col.isna().any():
                raise ValueError(
                    f"column {i} contains non-numeric or missing values; cannot convert to int"
                )
            arrays.append(col.astype(np.int64).to_numpy())
        else:
            arrays.append(col.astype(np.float64).to_numpy())
    return arrays


class BenchmarkResults:
    def __init__(self, arrays):
        (
            self.n,
            self.iterations,
            self.cpuTimes,
            self.cpuVtkmTimes,
            self.gpuVtkmTimes,
            self.gpuVtkmSpeedup,
            self.algoSpeedup,
            self.gpuSpeedup,
        ) = arrays

    @classmethod
    def from_sheet(cls, pathToSheet, startRow, endRow, colRange):
        arrays = readInResults(pathToSheet, startRow, endRow, colRange)
        return cls(arrays)

    def print(self):
        for name, attr in (
            ("n", "n"),
            ("iterations", "iterations"),
            ("cpuTimes", "cpuTimes"),
            ("cpuVtkmTimes", "cpuVtkmTimes"),
            ("gpuVtkmTimes", "gpuVtkmTimes"),
            ("gpuVtkmSpeedup", "gpuVtkmSpeedup"),
            ("algoSpeedup", "algoSpeedup"),
            ("gpuSpeedup", "gpuSpeedup"),
        ):
            print(f"{name}:")
            print(getattr(self, attr))


def plot_benchmarks(isoFull, isoPartial, clipFull, clipPartial):

    fig, axs = plt.subplots(1, 2, figsize=(14, 5), sharex="col")

    TITLE_FONT = {"family": "DejaVu Sans", "size": 18, "weight": "bold"}
    LABEL_FONT = {"family": "DejaVu Sans", "size": 16}
    TICK_FONT = {"family": "DejaVu Sans", "size": 14}

    labels = [
        "IsoSurface (original),  CPU",
        "IsoSurface (Viskores), CPU",
        "IsoSurface (Viskores), GPU",
        "Clip (original),  CPU",
        "Clip (Viskores), CPU",
        "Clip (Viskores), GPU",
    ]
    colors = ["#002772", "#009e24", "#FFBC00"]
    linesStyleFull = "-"
    linesStylePartial = "--"
    linewidth = 4

    suffix = " (rerun)"

    isoX = isoFull.n * isoFull.n * isoFull.n
    clipX = clipFull.n

    # PLOT LEFT
    axs[0].plot(
        isoX,
        isoFull.cpuTimes,
        linestyle=linesStyleFull,
        color=colors[0],
        label=labels[0],
        linewidth=linewidth,
    )
    """ axs[0].plot(
        isoX,
        isoPartial.cpuTimes,
        linestyle=linesStylePartial,
        color=colors[0],
        label=labels[0] + suffix,
        linewidth=linewidth,
    ) """

    axs[0].plot(
        isoX,
        isoFull.cpuVtkmTimes,
        linestyle=linesStyleFull,
        color=colors[1],
        label=labels[1],
        linewidth=linewidth,
    )
    """ axs[0].plot(
        isoX,
        isoPartial.cpuVtkmTimes,
        linestyle=linesStylePartial,
        color=colors[1],
        label=labels[1] + suffix,
        linewidth=linewidth,
    ) """

    axs[0].plot(
        isoX,
        isoFull.gpuVtkmTimes,
        linestyle=linesStyleFull,
        color=colors[2],
        label=labels[2],
        linewidth=linewidth,
    )
    """ axs[0].plot(
        isoX,
        isoPartial.gpuVtkmTimes,
        linestyle=linesStylePartial,
        color=colors[2],
        label=labels[2] + suffix,
        linewidth=linewidth,
    ) """

    # PLOT RIGHT
    axs[1].plot(
        clipX,
        clipFull.cpuTimes,
        linestyle=linesStyleFull,
        color=colors[0],
        label=labels[3],
        linewidth=linewidth,
    )
    """ axs[1].plot(
        clipX,
        clipPartial.cpuTimes,
        linestyle=linesStylePartial,
        color=colors[0],
        label=labels[3] + suffix,
        linewidth=linewidth,
    ) """

    axs[1].plot(
        clipX,
        clipFull.cpuVtkmTimes,
        linestyle=linesStyleFull,
        color=colors[1],
        label=labels[4],
        linewidth=linewidth,
    )
    """ axs[1].plot(
        clipX,
        clipPartial.cpuVtkmTimes,
        linestyle=linesStylePartial,
        color=colors[1],
        label=labels[4] + suffix,
        linewidth=linewidth,
    ) """

    axs[1].plot(
        clipX,
        clipFull.gpuVtkmTimes,
        linestyle=linesStyleFull,
        color=colors[2],
        label=labels[5],
        linewidth=linewidth,
    )
    """ axs[1].plot(
        clipX,
        clipPartial.gpuVtkmTimes,
        linestyle=linesStylePartial,
        color=colors[2],
        label=labels[5] + suffix,
        linewidth=linewidth,
    )
    """
    # formatting
    axs[0].set_title("IsoSurface Benchmark Runtimes", fontdict=TITLE_FONT)
    axs[1].set_title("Clip Benchmark Runtimes", fontdict=TITLE_FONT)

    for i in range(2):
        axs[i].set_xscale("log", base=10)  # linthresh controls linear region around 0

        # nicer ticks/formatting (base-10 major ticks + minor ticks)
        axs[i].xaxis.set_major_locator(LogLocator(base=10))
        axs[i].xaxis.set_minor_locator(LogLocator(base=10, subs=(2,3,4,5,6,7,8,9)))
        axs[i].xaxis.set_major_formatter(LogFormatter(base=10))
        axs[i].tick_params(which='major', length=6)
        axs[i].tick_params(which='minor', length=4) 
        

        axs[i].set_xlim(max(isoX[0], clipX[0]), max(isoX[-1], clipX[-1]))
        axs[i].set_ylim(0, 14)
        axs[i].set_xlabel("Number of Cells in Grid", fontdict=LABEL_FONT)
        axs[i].set_ylabel("Time [s]", fontdict=LABEL_FONT)
        axs[i].grid(True)

        axs[i].legend(loc="best", fontsize=TICK_FONT["size"])

        axs[i].tick_params(axis="both", which="major", labelsize=TICK_FONT["size"])
        for lbl in axs[i].get_xticklabels() + axs[i].get_yticklabels():
            lbl.set_fontfamily(TICK_FONT["family"])

    plt.tight_layout()
    outname = "benchmark_runtimes.svg"
    fig.savefig(outname, format="svg")
    print(f"Saved figure to {outname}")
    #plt.show()


if __name__ == "__main__":
    pathToSheet = "/home/hpcsmalh/Benchmarking/iso/overview.ods"

    isoFull = BenchmarkResults.from_sheet(pathToSheet, 5, 13, "A:H")
    isoPartial = BenchmarkResults.from_sheet(pathToSheet, 22, 30, "A:H")
    clipFull = BenchmarkResults.from_sheet(pathToSheet, 5, 18, "J:Q")
    clipPartial = BenchmarkResults.from_sheet(pathToSheet, 22, 35, "J:Q")

    plot_benchmarks(isoFull, isoPartial, clipFull, clipPartial)
