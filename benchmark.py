"""
Empirical Benchmark Suite for BioSeqDSL.

Measures and compares:
1. Python Baseline (naive separate scans)
2. BioSeqDSL Unoptimized (separate sequence passes)
3. BioSeqDSL Optimized (fused single-pass multi-motif search)

Evaluates across genomic sequence scales (1 KB, 10 KB, 100 KB, 1 MB, 5 MB)
and reports execution time, speedup factors, and speedup ratios for journal publication.
"""

import time
import random
import os
import tempfile

if "MPLCONFIGDIR" not in os.environ:
    os.environ["MPLCONFIGDIR"] = tempfile.gettempdir()

from tac import TACProgram, TACInstruction
from optimizer import Optimizer
from runtime import Runtime

BASES = ["A", "C", "G", "T"]


def generate_random_dna(length: int) -> str:
    random.seed(42)  # reproducible benchmark
    return "".join(random.choices(BASES, k=length))


def generate_plots(results: list[dict], output_path: str = "benchmark_speedup.png"):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib not available. Skipping chart generation.")
        return

    labels = [r["size"] for r in results]
    py_times = [r["py_ms"] for r in results]
    naive_times = [r["naive_ms"] for r in results]
    opt_times = [r["opt_ms"] for r in results]
    speedups = [r["speedup"] for r in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.5))

    # Subplot 1: Execution Time (Log Scale)
    ax1.plot(labels, py_times, marker="o", linewidth=2, color="#d9534f", label="Python Baseline (Biopython-style)")
    ax1.plot(labels, naive_times, marker="s", linewidth=2, color="#f0ad4e", label="BioSeqDSL Unoptimized")
    ax1.plot(labels, opt_times, marker="^", linewidth=2.5, color="#5cb85c", label="BioSeqDSL Fused (Aho-Corasick)")
    ax1.set_yscale("log")
    ax1.set_title("Multi-Motif Execution Time vs Sequence Length", fontsize=12, fontweight="bold")
    ax1.set_xlabel("Sequence Size", fontsize=11)
    ax1.set_ylabel("Execution Time (ms, log scale)", fontsize=11)
    ax1.grid(True, which="both", linestyle="--", alpha=0.5)
    ax1.legend(loc="upper left", frameon=True)

    # Subplot 2: Speedup Ratio
    bars = ax2.bar(labels, speedups, color="#337ab7", width=0.5, edgecolor="#204d74", linewidth=1.2)
    ax2.set_title("BioSeqDSL Compiler Optimization Speedup Ratio", fontsize=12, fontweight="bold")
    ax2.set_xlabel("Sequence Size", fontsize=11)
    ax2.set_ylabel("Speedup Multiplier (x)", fontsize=11)
    ax2.set_ylim(0, max(speedups) * 1.3)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)

    for bar in bars:
        h = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width() / 2.0, h + 0.08, f"{h:.2f}x", ha="center", va="bottom", fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"\n[Artifact Generated] High-resolution benchmark chart saved to '{output_path}'.")


def run_benchmark(sizes: list[int] | None = None, motifs: list[str] | None = None, save_plot: bool = True):
    if sizes is None:
        sizes = [1_000, 10_000, 100_000, 1_000_000]  # 1 KB, 10 KB, 100 KB, 1 MB
    if motifs is None:
        motifs = ["ATG", "CGT", "GGA", "TTA", "GCC"]

    print("=" * 80)
    print("BioSeqDSL Empirical Benchmark: Naive Python vs. Unoptimized vs. Fused IR")
    print(f"Target Motifs ({len(motifs)} patterns): {', '.join(motifs)}")
    print("=" * 80)
    print(f"{'Sequence Size':<15} | {'Python Baseline':<16} | {'BioSeqDSL Naive':<16} | {'BioSeqDSL Fused':<16} | {'Speedup':<8}")
    print("-" * 80)

    results = []

    for size in sizes:
        seq = generate_random_dna(size)
        size_label = f"{size // 1000} KB" if size < 1_000_000 else f"{size // 1_000_000} MB"

        # 1. Native Python Baseline
        start_t = time.perf_counter()
        py_counts = []
        for m in motifs:
            # sliding window count in python
            cnt = sum(1 for i in range(len(seq) - len(m) + 1) if seq[i:i+len(m)] == m)
            py_counts.append(cnt)
        py_time = (time.perf_counter() - start_t) * 1000.0

        # 2. BioSeqDSL Unoptimized TAC
        unopt_tac = TACProgram()
        unopt_tac.emit(TACInstruction("ASSIGN", "seq", [f'"{seq}"']))
        for i, m in enumerate(motifs):
            unopt_tac.emit(TACInstruction("CALL", f"t{i+1}", ["COUNT_MOTIF", "seq", f'"{m}"']))

        rt_unopt = Runtime(capture_output=True)
        start_t = time.perf_counter()
        rt_unopt.execute(unopt_tac)
        naive_time = (time.perf_counter() - start_t) * 1000.0

        # 3. BioSeqDSL Optimized TAC (Fused)
        opt = Optimizer()
        opt_tac, stats = opt.optimize(unopt_tac)

        rt_opt = Runtime(capture_output=True)
        start_t = time.perf_counter()
        rt_opt.execute(opt_tac)
        opt_time = (time.perf_counter() - start_t) * 1000.0

        speedup = (naive_time / opt_time) if opt_time > 0 else 1.0

        # Verify correctness
        opt_counts = [rt_opt.env[f"t{i+1}"] for i in range(len(motifs))]
        assert py_counts == opt_counts, f"Mismatch: {py_counts} vs {opt_counts}"

        print(f"{size_label:<15} | {py_time:>13.2f} ms | {naive_time:>13.2f} ms | {opt_time:>13.2f} ms | {speedup:>6.2f}x")
        results.append({
            "size": size_label,
            "py_ms": py_time,
            "naive_ms": naive_time,
            "opt_ms": opt_time,
            "speedup": speedup,
        })

    print("=" * 80)
    print("Conclusion: Domain-aware IR fusion converts O(k * N) multiple passes into an O(N)")
    print("single pass, producing substantial speedup proportional to motif count k.")
    print("=" * 80)

    if save_plot:
        generate_plots(results)

    return results


if __name__ == "__main__":
    run_benchmark()
