"""
BioSeqDSL Compiler & Execution Pipeline CLI.

Complete end-to-end compiler for DNA/RNA sequence analysis:
1. Lexical Analysis (Tokens)
2. Syntax Analysis (AST)
3. Semantic Analysis & Biological Type Checking (Symbol Table & Types)
4. Intermediate Code Generation (Three-Address Code / TAC)
5. Bioinformatics-Aware Optimization (Fused Motif Search)
6. Target Code Execution / Runtime

Usage:
    python main.py examples/01_length.bioseq
    python main.py examples/08_multiple_motif_searches.bioseq --all
    python main.py --benchmark
"""

import sys
import argparse

from lexer import Lexer
from parser import Parser
from ast_nodes import pretty
from semantic_analyzer import SemanticAnalyzer
from tac import TACGenerator
from optimizer import Optimizer
from runtime import Runtime
from errors import CompilerError


def run_pipeline(
    path: str,
    show_tokens: bool = False,
    show_ast: bool = False,
    show_symbols: bool = False,
    show_tac: bool = False,
    show_opt: bool = False,
    no_opt: bool = False,
    check_only: bool = False,
) -> int:
    try:
        with open(path, "r") as f:
            source = f.read()
    except OSError as err:
        print(f"Error reading file '{path}': {err}", file=sys.stderr)
        return 1

    source_lines = source.splitlines()

    try:
        # Phase 1: Lexical Analysis
        tokens = Lexer(source).tokenize()
        if show_tokens:
            print("=== 1. TOKENS ===")
            for tok in tokens:
                print(f"  {tok}")
            print()

        # Phase 2: Syntax Analysis
        program = Parser(tokens).parse_program()
        if show_ast:
            print("=== 2. ABSTRACT SYNTAX TREE (AST) ===")
            print(pretty(program))
            print()

        # Phase 3: Semantic Analysis & Biological Type Checking
        analyzer = SemanticAnalyzer()
        symbol_table = analyzer.analyze(program)
        if show_symbols:
            print("=== 3. SYMBOL TABLE & TYPES ===")
            print(symbol_table)
            print()

        # Phase 4: Intermediate Code Generation (TAC)
        tac_prog = TACGenerator().generate(program)
        if show_tac:
            print("=== 4. THREE-ADDRESS CODE (UNOPTIMIZED TAC) ===")
            print(tac_prog)
            print()

        # Phase 5: Domain-Specific Optimization
        if not no_opt:
            optimizer = Optimizer()
            opt_tac, stats = optimizer.optimize(tac_prog)
            target_tac = opt_tac
            if show_opt:
                print("=== 5. OPTIMIZED TAC (FUSED MOTIF SEARCH) ===")
                print(target_tac)
                print(f"\n{stats}\n")
        else:
            target_tac = tac_prog
            if show_opt:
                print("=== 5. OPTIMIZATION DISABLED (--no-opt) ===\n")

        # Phase 6: Code Execution / Runtime
        if not check_only:
            if show_tokens or show_ast or show_symbols or show_tac or show_opt:
                print("=== 6. PROGRAM OUTPUT ===")
            runtime = Runtime(capture_output=False)
            runtime.execute(target_tac)
        else:
            print(f"Program is valid. ({path})")

        return 0

    except CompilerError as e:
        print(e.render(source_lines), file=sys.stderr)
        return 1


def main():
    arg_parser = argparse.ArgumentParser(
        description="BioSeqDSL Compiler & Execution Environment for DNA/RNA Bioinformatics"
    )
    arg_parser.add_argument("file", nargs="?", help="Path to a .bioseq source file")
    arg_parser.add_argument("--tokens", action="store_true", help="Print the token stream")
    arg_parser.add_argument("--ast", action="store_true", help="Print the AST")
    arg_parser.add_argument("--symbols", action="store_true", help="Print the symbol table")
    arg_parser.add_argument("--tac", action="store_true", help="Print unoptimized Three-Address Code")
    arg_parser.add_argument("--opt", action="store_true", help="Print optimized Three-Address Code")
    arg_parser.add_argument("--no-opt", action="store_true", help="Disable domain-specific optimizations")
    arg_parser.add_argument("--check-only", action="store_true", help="Validate syntax and types without running")
    arg_parser.add_argument("--all", action="store_true", help="Display all intermediate stages and output")
    arg_parser.add_argument("--benchmark", action="store_true", help="Run the empirical performance benchmark suite")

    args = arg_parser.parse_args()

    if args.benchmark:
        from benchmark import run_benchmark
        run_benchmark()
        return

    if not args.file:
        arg_parser.print_help()
        sys.exit(1)

    show_tokens = args.tokens or args.all
    show_ast = args.ast or args.all
    show_symbols = args.symbols or args.all
    show_tac = args.tac or args.all
    show_opt = args.opt or args.all

    sys.exit(run_pipeline(
        args.file,
        show_tokens=show_tokens,
        show_ast=show_ast,
        show_symbols=show_symbols,
        show_tac=show_tac,
        show_opt=show_opt,
        no_opt=args.no_opt,
        check_only=args.check_only,
    ))


if __name__ == "__main__":
    main()
