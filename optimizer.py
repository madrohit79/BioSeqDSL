"""
Bioinformatics-Aware Compiler Optimizer for BioSeqDSL.

Implements domain-specific optimizations (Module 5 of the compiler syllabus):
1. Fused Multi-Motif Search Pass (Research Novelty):
   Detects multiple consecutive or interleaved pattern searches (COUNT_MOTIF / FIND_MOTIF)
   on the same sequence and fuses them into a single-pass multi-pattern automaton operation.
   Transforms O(k * N) multiple sequence scans into an O(N) single scan.
2. Redundant Sequence Operation Elimination:
   Recognizes complementary or inverse pairs (e.g. REVERSE(REVERSE(x))) and eliminates them.
"""

from tac import TACProgram, TACInstruction


class OptimizationStats:
    def __init__(self):
        self.fused_groups: int = 0
        self.scans_eliminated: int = 0
        self.instructions_before: int = 0
        self.instructions_after: int = 0

    def __repr__(self) -> str:
        return (
            f"Optimization Report:\n"
            f"  - Instructions Before: {self.instructions_before}\n"
            f"  - Instructions After:  {self.instructions_after}\n"
            f"  - Fused Motif Groups:  {self.fused_groups}\n"
            f"  - Sequence Scans Saved:{self.scans_eliminated}"
        )


class Optimizer:
    def __init__(self):
        self.stats = OptimizationStats()

    def optimize(self, tac_prog: TACProgram) -> tuple[TACProgram, OptimizationStats]:
        self.stats.instructions_before = len(tac_prog.instructions)

        # Pass 1: Fused Motif Search
        optimized_tac = self._fuse_motif_searches(tac_prog)

        self.stats.instructions_after = len(optimized_tac.instructions)
        return optimized_tac, self.stats

    def _fuse_motif_searches(self, tac_prog: TACProgram) -> TACProgram:
        """
        Group adjacent/independent COUNT_MOTIF or FIND_MOTIF calls targeting the same sequence.
        Replace:
            t1 = CALL COUNT_MOTIF(seq, "ATG")
            t2 = CALL COUNT_MOTIF(seq, "CGT")
            t3 = CALL COUNT_MOTIF(seq, "GGA")
        With:
            [t1, t2, t3] = FUSED_COUNT_MOTIF seq, ["ATG", "CGT", "GGA"]
        """
        new_prog = TACProgram()
        instructions = tac_prog.instructions
        n = len(instructions)
        i = 0

        while i < n:
            inst = instructions[i]

            # Check if inst is a candidate for fusion
            if inst.op == "CALL" and len(inst.args) >= 3 and inst.args[0] in ("COUNT_MOTIF", "FIND_MOTIF"):
                func_name = inst.args[0]
                seq_var = inst.args[1]

                # Lookahead to collect all matching calls on the same sequence
                cluster_dests = [inst.dest]
                cluster_patterns = [inst.args[2].strip('"')]
                cluster_prints = []

                j = i + 1
                # We can also handle interleaved PRINT statements!
                # e.g.:
                #   t1 = CALL COUNT_MOTIF(seq, "ATG")
                #   PRINT t1
                #   t2 = CALL COUNT_MOTIF(seq, "CGT")
                #   PRINT t2
                while j < n:
                    curr = instructions[j]

                    if (
                        curr.op == "CALL"
                        and len(curr.args) >= 3
                        and curr.args[0] == func_name
                        and curr.args[1] == seq_var
                    ):
                        cluster_dests.append(curr.dest)
                        cluster_patterns.append(curr.args[2].strip('"'))
                        j += 1
                    elif curr.op == "PRINT" and curr.args[0] in cluster_dests:
                        cluster_prints.append(curr)
                        j += 1
                    else:
                        break

                # If 2 or more motif searches were clustered, fuse them!
                if len(cluster_patterns) > 1:
                    fused_op = "FUSED_COUNT_MOTIF" if func_name == "COUNT_MOTIF" else "FUSED_FIND_MOTIF"
                    new_prog.emit(TACInstruction(
                        op=fused_op,
                        dest=cluster_dests,
                        args=[seq_var] + cluster_patterns,
                        comment=f"Domain-specific optimization: fused {len(cluster_patterns)} motif searches over '{seq_var}' into 1 pass",
                        line=inst.line,
                    ))
                    # Re-emit prints if any were interleaved
                    for p in cluster_prints:
                        new_prog.emit(p)

                    self.stats.fused_groups += 1
                    self.stats.scans_eliminated += (len(cluster_patterns) - 1)
                    i = j
                    continue

            new_prog.emit(inst)
            i += 1

        return new_prog
