"""
Runtime Execution Engine for BioSeqDSL.

Executes Three-Address Code (TAC) representations:
- Supports standard and degenerate sequence biology operations (transcription, translation, GC content, etc.)
- Full IUPAC nucleotide degenerate alphabet support
- Complete Aho-Corasick Multi-Pattern Automaton for fused motif scanning (O(N))
- Control flow branching (IF / ELSE, WHILE loops via program counter jumps)
- FASTA file I/O (LOAD_FASTA, SAVE_FASTA)
- Formatted biological output rendering.
"""

from typing import Any
from collections import deque
import os
from tac import TACProgram, TACInstruction

CODON_TABLE = {
    # Phenylalanine
    "UUU": "F", "UUC": "F",
    # Leucine
    "UUA": "L", "UUG": "L", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    # Isoleucine
    "AUU": "I", "AUC": "I", "AUA": "I",
    # Methionine (Start)
    "AUG": "M",
    # Valine
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    # Serine
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S", "AGU": "S", "AGC": "S",
    # Proline
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Threonine
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Alanine
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyrosine
    "UAU": "Y", "UAC": "Y",
    # Stop Codons
    "UAA": "*", "UAG": "*", "UGA": "*",
    # Histidine
    "CAU": "H", "CAC": "H",
    # Glutamine
    "CAA": "Q", "CAG": "Q",
    # Asparagine
    "AAU": "N", "AAC": "N",
    # Lysine
    "AAA": "K", "AAG": "K",
    # Aspartate
    "GAU": "D", "GAC": "D",
    # Glutamate
    "GAA": "E", "GAG": "E",
    # Cysteine
    "UGU": "C", "UGC": "C",
    # Tryptophan
    "UGG": "W",
    # Arginine
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Glycine
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Complete IUPAC Degenerate Nucleotide Complements:
# A <-> T (or U), C <-> G
# R (A/G) <-> Y (T/C)
# S (C/G) <-> S (G/C), W (A/T) <-> W (T/A)
# K (G/T) <-> M (A/C)
# B (C/G/T) <-> V (A/C/G), D (A/G/T) <-> H (A/C/T)
# N <-> N
DNA_IUPAC_FROM = "ACGTURYWSKMBDHVNacgturywskmbdhvn"
DNA_IUPAC_TO   = "TGCAAYRWSMKVHDBNtgcaayrwsmkvhdbn"
DNA_COMPLEMENT = str.maketrans(DNA_IUPAC_FROM, DNA_IUPAC_TO)

RNA_IUPAC_FROM = "ACGTURYWSKMBDHVNacgturywskmbdhvn"
RNA_IUPAC_TO   = "UGCAAYRWSMKVHDBNugcaayrwsmkvhdbn"
RNA_COMPLEMENT = str.maketrans(RNA_IUPAC_FROM, RNA_IUPAC_TO)


class AhoCorasickNode:
    def __init__(self):
        self.children: dict[str, AhoCorasickNode] = {}
        self.fail: AhoCorasickNode | None = None
        self.patterns: list[str] = []


class AhoCorasick:
    """Aho-Corasick Multi-Pattern Automaton for O(N + sum(M_i)) search."""
    def __init__(self, patterns: list[str]):
        self.root = AhoCorasickNode()
        self.patterns = patterns
        self._build_trie()
        self._build_failure_links()

    def _build_trie(self):
        for pat in self.patterns:
            node = self.root
            for ch in pat:
                if ch not in node.children:
                    node.children[ch] = AhoCorasickNode()
                node = node.children[ch]
            node.patterns.append(pat)

    def _build_failure_links(self):
        queue = deque()
        for ch, child in self.root.children.items():
            child.fail = self.root
            queue.append(child)

        while queue:
            curr = queue.popleft()
            for ch, child in curr.children.items():
                fail_node = curr.fail
                while fail_node is not None and ch not in fail_node.children:
                    fail_node = fail_node.fail
                child.fail = fail_node.children[ch] if fail_node else self.root
                if child.fail.patterns:
                    child.patterns.extend(child.fail.patterns)
                queue.append(child)

    def count(self, text: str) -> dict[str, int]:
        counts = {p: 0 for p in self.patterns}
        node = self.root
        for ch in text:
            while node is not None and ch not in node.children:
                node = node.fail
            if node is None:
                node = self.root
                continue
            node = node.children[ch]
            for p in node.patterns:
                counts[p] += 1
        return counts

    def search(self, text: str) -> dict[str, list[int]]:
        hits: dict[str, list[int]] = {p: [] for p in self.patterns}
        node = self.root
        for idx, ch in enumerate(text):
            while node is not None and ch not in node.children:
                node = node.fail
            if node is None:
                node = self.root
                continue
            node = node.children[ch]
            for p in node.patterns:
                start_pos = (idx - len(p) + 1) + 1  # 1-based index
                hits[p].append(start_pos)
        return hits


class Runtime:
    def __init__(self, capture_output: bool = False):
        self.env: dict[str, Any] = {}
        self.capture_output = capture_output
        self.output_lines: list[str] = []

    def execute(self, tac_prog: TACProgram) -> list[str]:
        self.env.clear()
        self.output_lines.clear()

        # Build label map for GOTO and IF_FALSE branching
        label_map: dict[str, int] = {}
        for idx, inst in enumerate(tac_prog.instructions):
            if inst.op == "LABEL":
                label_map[inst.args[0]] = idx

        instructions = tac_prog.instructions
        n = len(instructions)
        pc = 0

        while pc < n:
            inst = instructions[pc]
            op = inst.op

            if op == "LABEL":
                pc += 1

            elif op == "GOTO":
                target = inst.args[0]
                pc = label_map[target]

            elif op == "IF_FALSE":
                cond_val = self._resolve_val(inst.args[0])
                target = inst.args[1]
                is_false = (cond_val is False) or (cond_val == "FALSE") or (cond_val == 0) or (not cond_val)
                if is_false:
                    pc = label_map[target]
                else:
                    pc += 1

            elif op == "BINARY_OP":
                op_name = inst.args[0]
                left = self._resolve_val(inst.args[1])
                right = self._resolve_val(inst.args[2])
                self.env[inst.dest] = self._eval_binary_op(op_name, left, right)
                pc += 1

            elif op == "UNARY_OP":
                op_name = inst.args[0]
                val = self._resolve_val(inst.args[1])
                if op_name == "NOT":
                    is_true = (val is True) or (val == "TRUE")
                    self.env[inst.dest] = not is_true
                pc += 1

            elif op == "ASSIGN":
                val = self._resolve_val(inst.args[0])
                self.env[inst.dest] = val
                pc += 1

            elif op == "CALL":
                func = inst.args[0]
                evaluated_args = [self._resolve_val(a) for a in inst.args[1:]]
                res = self._call_builtin(func, evaluated_args)
                self.env[inst.dest] = res
                pc += 1

            elif op == "FUSED_COUNT_MOTIF":
                seq = self._resolve_val(inst.args[0])
                patterns = inst.args[1:]
                counts = self._fused_count_motif(seq, patterns)
                dests = inst.dest if isinstance(inst.dest, list) else [inst.dest]
                for dest, count in zip(dests, counts):
                    self.env[dest] = count
                pc += 1

            elif op == "FUSED_FIND_MOTIF":
                seq = self._resolve_val(inst.args[0])
                patterns = inst.args[1:]
                hits = self._fused_find_motif(seq, patterns)
                dests = inst.dest if isinstance(inst.dest, list) else [inst.dest]
                for dest, hit_list in zip(dests, hits):
                    self.env[dest] = hit_list
                pc += 1

            elif op == "PRINT":
                val = self._resolve_val(inst.args[0])
                out_str = self._format_output(val)
                self.output_lines.append(out_str)
                if not self.capture_output:
                    print(out_str)
                pc += 1

            else:
                pc += 1

        return self.output_lines

    def _resolve_val(self, operand: Any) -> Any:
        if isinstance(operand, (int, float, bool, list)):
            return operand
        operand = str(operand).strip()

        # Boolean literals
        if operand in ("TRUE", "true"):
            return True
        if operand in ("FALSE", "false"):
            return False

        # String literal: "..."
        if operand.startswith('"') and operand.endswith('"'):
            return operand[1:-1]

        # Integer literal
        if operand.isdigit() or (operand.startswith('-') and operand[1:].isdigit()):
            return int(operand)

        # Float literal
        try:
            if "." in operand:
                return float(operand)
        except ValueError:
            pass

        # Environment lookup
        if operand in self.env:
            return self.env[operand]

        return operand

    def _eval_binary_op(self, op: str, left: Any, right: Any) -> Any:
        if op == "+":
            return left + right
        if op == "-":
            return left - right
        if op == "<":
            return left < right
        if op == ">":
            return left > right
        if op == "<=":
            return left <= right
        if op == ">=":
            return left >= right
        if op == "==":
            return left == right
        if op == "!=":
            return left != right
        if op == "AND":
            b_left = (left is True) or (left == "TRUE")
            b_right = (right is True) or (right == "TRUE")
            return b_left and b_right
        if op == "OR":
            b_left = (left is True) or (left == "TRUE")
            b_right = (right is True) or (right == "TRUE")
            return b_left or b_right
        raise ValueError(f"Unknown binary operator: {op}")

    def _call_builtin(self, func: str, args: list[Any]) -> Any:
        seq = str(args[0])

        if func == "LENGTH":
            return len(seq)

        if func == "GC_CONTENT":
            if not seq:
                return 0.0
            seq_upper = seq.upper()
            gc = seq_upper.count("G") + seq_upper.count("C")
            return round((gc / len(seq)) * 100.0, 2)

        if func == "REVERSE":
            return seq[::-1]

        if func == "COMPLEMENT":
            trans = RNA_COMPLEMENT if "U" in seq.upper() else DNA_COMPLEMENT
            return seq.translate(trans)

        if func == "REVERSE_COMPLEMENT":
            trans = RNA_COMPLEMENT if "U" in seq.upper() else DNA_COMPLEMENT
            return seq.translate(trans)[::-1]

        if func == "TRANSCRIBE":
            return seq.replace("T", "U").replace("t", "u")

        if func == "TRANSLATE":
            rna = seq.upper()
            protein = []
            for i in range(0, len(rna) - 2, 3):
                codon = rna[i:i+3]
                aa = CODON_TABLE.get(codon, "?")
                if aa == "*":
                    break
                protein.append(aa)
            return "".join(protein)

        if func == "FIND_MOTIF":
            pattern = str(args[1]).upper()
            seq_upper = seq.upper()
            ac = AhoCorasick([pattern])
            hits = ac.search(seq_upper)
            return hits[pattern]

        if func == "COUNT_MOTIF":
            pattern = str(args[1]).upper()
            seq_upper = seq.upper()
            ac = AhoCorasick([pattern])
            counts = ac.count(seq_upper)
            return counts[pattern]

        if func == "LOAD_FASTA":
            filepath = str(args[0])
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"FASTA file not found: {filepath}")
            lines = []
            with open(filepath, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith(">") or line.startswith(";"):
                        continue
                    lines.append(line)
            return "".join(lines).upper()

        if func == "SAVE_FASTA":
            filepath = str(args[0])
            sequence = str(args[1])
            header = str(args[2])
            with open(filepath, "w") as f:
                f.write(f">{header}\n")
                for i in range(0, len(sequence), 70):
                    f.write(sequence[i:i+70] + "\n")
            return len(sequence)

        raise ValueError(f"Unknown builtin function: {func}")

    def _fused_count_motif(self, seq: str, patterns: list[str]) -> list[int]:
        """Aho-Corasick automaton single-pass multi-pattern search."""
        seq_upper = seq.upper()
        pats_upper = [p.upper() for p in patterns]
        ac = AhoCorasick(pats_upper)
        counts = ac.count(seq_upper)
        return [counts[p] for p in pats_upper]

    def _fused_find_motif(self, seq: str, patterns: list[str]) -> list[list[int]]:
        """Aho-Corasick automaton single-pass multi-pattern location collection."""
        seq_upper = seq.upper()
        pats_upper = [p.upper() for p in patterns]
        ac = AhoCorasick(pats_upper)
        hits = ac.search(seq_upper)
        return [hits[p] for p in pats_upper]

    def _format_output(self, val: Any) -> str:
        if isinstance(val, bool):
            return "TRUE" if val else "FALSE"
        if isinstance(val, float):
            return f"{val:.2f}%" if val <= 100.0 else f"{val:.2f}"
        if isinstance(val, list):
            return "[" + ", ".join(str(x) for x in val) + "]"
        return str(val)
