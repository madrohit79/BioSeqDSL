# BioSeqDSL: A Domain-Specific Language & Optimizing Compiler for Bioinformatics

> **Project Title:** Compiler for a Tiny Domain-Specific Language for Bioinformatics  
> **Team 29:** Rohit M (24BKT0055), Shivam Aditya R (24BKT0007), K Vaishnav S Reddy (24BKT0047)  
> **Target Outcome:** Journal Publication & Compiler Design Course Project

---

## 1. Project Overview & Motivation

General-purpose programming languages (Python, R, Java) and bioinformatics libraries (Biopython, BioPerl) treat biological sequences as generic strings. Consequently, elementary sequence operations (transcription, translation, GC content, motif searches) require developers to reconstruct biological logic using low-level string slicing and loops. This creates an abstraction gap and prevents compilers from performing domain-specific optimizations.

**BioSeqDSL** is a domain-specific language and optimizing compiler designed specifically for DNA/RNA sequence manipulation. It makes biological meaning a first-class compiler construct, enforcing biological type constraints ahead of execution and applying **Domain-Aware Intermediate Representation (IR) optimizations** such as **Fused Multi-Motif Search**.

---

## 2. Compiler Architecture & Syllabus Mapping

BioSeqDSL implements the entire 7-module compiler pipeline:

```
  Source Code (.bioseq)
           │
           ▼
  ┌─────────────────┐
  │ 1. Lexer        │  Module 1: Lexical Analysis (Tokens, Line/Col spans)
  └────────┬────────┘
           │ Tokens
           ▼
  ┌─────────────────┐
  │ 2. Parser       │  Module 2: LL(1) Recursive Descent Parser (AST)
  └────────┬────────┘
           │ AST
           ▼
  ┌─────────────────┐
  │ 3. Semantic     │  Module 3: Type Checking & Alphabet Validation
  │    Analyzer     │            (DNA/RNA alphabets, TRANSCRIBE/TRANSLATE constraints)
  └────────┬────────┘
           │ Checked AST + Symbol Table
           ▼
  ┌─────────────────┐
  │ 4. TAC Gen      │  Module 4: Intermediate Code Generation (Three-Address Code)
  └────────┬────────┘
           │ Unoptimized TAC
           ▼
  ┌─────────────────┐
  │ 5. Optimizer    │  Module 5: Domain-Aware Optimization (Fused Motif Automaton)
  └────────┬────────┘
           │ Optimized TAC
           ▼
  ┌─────────────────┐
  │ 6. Runtime / VM │  Module 6 & 7: Code Execution & DSL Evaluation
  └────────┬────────┘
           │
           ▼
     Biological Output
```

---

## 3. The Research Novelty: Fused Multi-Motif Search

When searching for multiple motifs over the same sequence, naive compilers and libraries perform $k$ separate scans ($O(k \cdot N)$):

```
PRINT COUNT_MOTIF(seq, "ATG");
PRINT COUNT_MOTIF(seq, "CGT");
PRINT COUNT_MOTIF(seq, "GGA");
```

### Unoptimized IR (Naive):
```text
002:  t1 = CALL COUNT_MOTIF(seq, "ATG")
003:  PRINT t1
004:  t2 = CALL COUNT_MOTIF(seq, "CGT")
005:  PRINT t2
006:  t3 = CALL COUNT_MOTIF(seq, "GGA")
007:  PRINT t3
```

### BioSeqDSL Optimized IR (Fused Pass):
The optimizer recognizes adjacent motif queries on the same sequence and fuses them into a single-pass prefix-dispatched search:
```text
002:  [t1, t2, t3] = FUSED_COUNT_MOTIF seq, ['ATG', 'CGT', 'GGA']
003:  PRINT t1
004:  PRINT t2
005:  PRINT t3
```

### Empirical Results (from `python main.py --benchmark`):
| Sequence Size | Python Baseline | BioSeqDSL Naive | BioSeqDSL Fused (Aho-Corasick) | Speedup Ratio |
| :--- | :--- | :--- | :--- | :--- |
| **1 KB** | 0.33 ms | 0.31 ms | 0.10 ms | **3.24x** |
| **10 KB** | 3.43 ms | 2.53 ms | 0.80 ms | **3.18x** |
| **100 KB** | 35.52 ms | 25.52 ms | 8.94 ms | **2.86x** |
| **1 MB** | 352.07 ms | 258.87 ms | 88.47 ms | **2.93x** |

Running the benchmark also generates a high-resolution publication figure `benchmark_speedup.png`.

---

## 4. Advanced Features (Review 3 Deliverables)

### 1. Control Flow Structures (`IF` / `ELSE`, `WHILE`)
BioSeqDSL supports conditional branching and loops with comparison (`<`, `>`, `<=`, `>=`, `==`, `!=`) and logical operators (`AND`, `OR`, `NOT`):

```bioseq
DNA seq = "GCGCGCGCGCGCGCGC";
FLOAT gc = GC_CONTENT(seq);

IF (gc > 60.0) {
    PRINT "High GC content sequence";
    PRINT gc;
} ELSE {
    PRINT "Standard GC sequence";
}

INTEGER countdown = 3;
WHILE (countdown > 0) {
    PRINT countdown;
    countdown = countdown - 1;
}
```

### 2. FASTA File Ingestion & Export (`LOAD_FASTA`, `SAVE_FASTA`)
Directly read and write standard NCBI/Ensembl multi-line FASTA files:

```bioseq
DNA genome = LOAD_FASTA("examples/sample_gene.fasta");
PRINT LENGTH(genome);
PRINT GC_CONTENT(genome);

DNA rcomp = REVERSE_COMPLEMENT(genome);
INTEGER written = SAVE_FASTA("output.fasta", rcomp, "Reverse Complement");
```

### 3. IUPAC Degenerate Nucleotide Alphabets
Supports ambiguous sequencing reads using IUPAC degenerate codes (`R` = A/G, `Y` = C/T, `S` = G/C, `W` = A/T, `K` = G/T, `M` = A/C, `B` = C/G/T, `D` = A/G/T, `H` = A/C/T, `V` = A/C/G, `N` = any nucleotide) with accurate biological complement transformations:

```bioseq
DNA degenerate_probe = "ATGNNNCRYWSKMBVHD";
PRINT COMPLEMENT(degenerate_probe);
// Output: TACNNNGYRWSMKVBDH
```

### 4. Aho-Corasick Multi-Pattern Automaton
Fused pattern matching is powered by a textbook Aho-Corasick finite state automaton with trie prefix links and failure transitions, guaranteeing optimal $O(N + \sum M_i)$ runtime across hundreds of simultaneous motifs.

---

## 5. Quick Start & CLI Usage

Run any BioSeqDSL program with:

```bash
# Execute program and print biological outputs
python main.py examples/05_transcribe_translate.bioseq

# Inspect all 6 compiler stages (Tokens, AST, Symbols, TAC, Opt TAC, Output)
python main.py examples/08_multiple_motif_searches.bioseq --all

# Run without optimizations
python main.py examples/08_multiple_motif_searches.bioseq --no-opt

# Run empirical benchmark suite
python main.py --benchmark
```

### Running Unit Tests

```bash
python tests/test_pipeline.py
# or with pytest:
pytest
```

---

## 6. Rich Diagnostic System

BioSeqDSL features source-aware error reporting with exact line, column, and caret spans:

```text
Error: Invalid Argument Type (line 2, column 18)
  STRING protein = TRANSLATE(seq);
                   ^^^^^^^^^
  TRANSLATE expects RNA, received DNA
  Hint: TRANSLATE operates on RNA. Did you forget TRANSCRIBE(seq) first?
```
