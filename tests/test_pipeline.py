"""
Comprehensive test suite for the complete BioSeqDSL compiler pipeline.

Covers:
1. Lexer (Tokens, source tracking, comments)
2. Parser (AST construction, grammar enforcement)
3. Symbol Table (Scope, undeclared/redeclaration checks)
4. Semantic Analyzer (Type checking, sequence alphabets, biological operation signatures)
5. Intermediate Code Generator (Three-Address Code lowering)
6. Optimizer (Domain-specific Fused Motif Search pass)
7. Runtime Execution (Correct biological results for all operations)

Can be run with:
    pytest
or:
    python tests/test_pipeline.py
"""

import sys
import os
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from lexer import Lexer
from parser import Parser
from tokens import TokenType
from ast_nodes import Program, Declaration, PrintStatement, FunctionCall
from symbol_table import build_symbol_table
from semantic_analyzer import SemanticAnalyzer
from tac import TACGenerator
from optimizer import Optimizer
from runtime import Runtime
from errors import (
    UndeclaredIdentifierError, RedeclarationError, ParseError,
    LexError, TypeMismatchError, InvalidArgumentTypeError, InvalidAlphabetError,
)


class TestBioSeqDSLPipeline(unittest.TestCase):

    # ---------- 1. Lexer tests ----------

    def test_lexer_tokenizes_declaration(self):
        tokens = Lexer('DNA seq = "ATGC";').tokenize()
        types = [t.type for t in tokens]
        self.assertEqual(types, [
            TokenType.DNA,
            TokenType.IDENTIFIER,
            TokenType.ASSIGN,
            TokenType.STRING_LITERAL,
            TokenType.SEMICOLON,
            TokenType.EOF,
        ])

    def test_lexer_tracks_line_and_column(self):
        tokens = Lexer('DNA seq = "ATGC";\nPRINT seq;').tokenize()
        print_tok = next(t for t in tokens if t.type == TokenType.PRINT)
        self.assertEqual(print_tok.line, 2)
        self.assertEqual(print_tok.column, 1)

    def test_lexer_rejects_illegal_character(self):
        with self.assertRaises(LexError):
            Lexer("DNA seq = @;").tokenize()

    def test_lexer_skips_comments(self):
        tokens = Lexer("// a comment\nPRINT 1;").tokenize()
        types = [t.type for t in tokens]
        self.assertEqual(types, [TokenType.PRINT, TokenType.INTEGER_LITERAL, TokenType.SEMICOLON, TokenType.EOF])

    # ---------- 2. Parser tests ----------

    def test_parser_builds_declaration_node(self):
        tokens = Lexer('DNA seq = "ATGC";').tokenize()
        program = Parser(tokens).parse_program()
        self.assertIsInstance(program, Program)
        self.assertEqual(len(program.statements), 1)
        decl = program.statements[0]
        self.assertIsInstance(decl, Declaration)
        self.assertEqual(decl.var_type, "DNA")
        self.assertEqual(decl.name, "seq")

    def test_parser_builds_function_call_node(self):
        tokens = Lexer('PRINT GC_CONTENT(seq);').tokenize()
        program = Parser(tokens).parse_program()
        stmt = program.statements[0]
        self.assertIsInstance(stmt, PrintStatement)
        self.assertIsInstance(stmt.value, FunctionCall)
        self.assertEqual(stmt.value.name, "GC_CONTENT")
        self.assertEqual(len(stmt.value.arguments), 1)

    def test_parser_rejects_missing_semicolon(self):
        tokens = Lexer('DNA seq = "ATGC"\nPRINT seq;').tokenize()
        with self.assertRaises(ParseError):
            Parser(tokens).parse_program()

    # ---------- 3. Symbol table tests ----------

    def test_symbol_table_declares_identifier(self):
        tokens = Lexer('DNA seq = "ATGC";').tokenize()
        program = Parser(tokens).parse_program()
        table = build_symbol_table(program)
        info = table.lookup("seq", 1, 1)
        self.assertEqual(info.var_type, "DNA")

    def test_symbol_table_catches_undeclared_identifier(self):
        tokens = Lexer('PRINT LENGTH(other);').tokenize()
        program = Parser(tokens).parse_program()
        with self.assertRaises(UndeclaredIdentifierError):
            build_symbol_table(program)

    def test_symbol_table_catches_redeclaration(self):
        tokens = Lexer('DNA seq = "ATGC";\nDNA seq = "TTTT";').tokenize()
        program = Parser(tokens).parse_program()
        with self.assertRaises(RedeclarationError):
            build_symbol_table(program)

    # ---------- 4. Semantic Analyzer & Type Checking tests ----------

    def test_semantic_rejects_dna_with_invalid_alphabet(self):
        tokens = Lexer('DNA seq = "ATGCXTAG";').tokenize()
        program = Parser(tokens).parse_program()
        with self.assertRaises(InvalidAlphabetError):
            SemanticAnalyzer().analyze(program)

    def test_semantic_rejects_rna_with_thymine(self):
        tokens = Lexer('RNA seq = "AUGCT";').tokenize()
        program = Parser(tokens).parse_program()
        with self.assertRaises(InvalidAlphabetError):
            SemanticAnalyzer().analyze(program)

    def test_semantic_rejects_type_mismatch_assignment(self):
        source = 'DNA seq = "ATGC";\nRNA rna = "AUGC";\nDNA copy = rna;'
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        with self.assertRaises(TypeMismatchError):
            SemanticAnalyzer().analyze(program)

    def test_semantic_rejects_translate_on_dna(self):
        source = 'DNA seq = "ATGGCC";\nSTRING prot = TRANSLATE(seq);'
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        with self.assertRaises(InvalidArgumentTypeError):
            SemanticAnalyzer().analyze(program)

    # ---------- 5. TAC & Optimization tests ----------

    def test_tac_generation_and_fused_optimization(self):
        source = '''
        DNA seq = "ATGCGTATGCCGATG";
        PRINT COUNT_MOTIF(seq, "ATG");
        PRINT COUNT_MOTIF(seq, "CGT");
        PRINT COUNT_MOTIF(seq, "GGA");
        '''
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        SemanticAnalyzer().analyze(program)
        tac = TACGenerator().generate(program)
        self.assertEqual(len(tac.instructions), 7)

        opt_tac, stats = Optimizer().optimize(tac)
        self.assertEqual(stats.fused_groups, 1)
        self.assertEqual(stats.scans_eliminated, 2)
        fused_inst = next(i for i in opt_tac.instructions if i.op == "FUSED_COUNT_MOTIF")
        self.assertEqual(fused_inst.args, ["seq", "ATG", "CGT", "GGA"])

    # ---------- 6. Runtime Execution tests ----------

    def test_runtime_transcribe_and_translate(self):
        source = '''
        DNA gene = "ATGGCCATTGTAATG";
        RNA transcript = TRANSCRIBE(gene);
        STRING protein = TRANSLATE(transcript);
        PRINT protein;
        '''
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        SemanticAnalyzer().analyze(program)
        tac = TACGenerator().generate(program)
        rt = Runtime(capture_output=True)
        lines = rt.execute(tac)
        self.assertEqual(lines, ["MAIVM"])

    def test_runtime_gc_content_and_reverse_complement(self):
        source = '''
        DNA seq = "ATGC";
        PRINT GC_CONTENT(seq);
        PRINT REVERSE_COMPLEMENT(seq);
        '''
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        SemanticAnalyzer().analyze(program)
        tac = TACGenerator().generate(program)
        rt = Runtime(capture_output=True)
        lines = rt.execute(tac)
        self.assertEqual(lines, ["50.00%", "GCAT"])

    # ---------- 7. Control Flow & Advanced Features (Review 3) ----------

    def test_control_flow_if_else(self):
        source = '''
        DNA rich = "GCGCGC";
        FLOAT gc = GC_CONTENT(rich);
        IF (gc > 50.0) {
            PRINT "high";
        } ELSE {
            PRINT "low";
        }
        '''
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        SemanticAnalyzer().analyze(program)
        tac = TACGenerator().generate(program)
        rt = Runtime(capture_output=True)
        lines = rt.execute(tac)
        self.assertEqual(lines, ["high"])

    def test_control_flow_while_loop(self):
        source = '''
        INTEGER i = 3;
        WHILE (i > 0) {
            PRINT i;
            i = i - 1;
        }
        '''
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        SemanticAnalyzer().analyze(program)
        tac = TACGenerator().generate(program)
        rt = Runtime(capture_output=True)
        lines = rt.execute(tac)
        self.assertEqual(lines, ["3", "2", "1"])

    def test_semantic_rejects_non_boolean_if_condition(self):
        source = '''
        DNA seq = "ATGC";
        IF (seq) {
            PRINT 1;
        }
        '''
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        with self.assertRaises(TypeMismatchError):
            SemanticAnalyzer().analyze(program)

    def test_iupac_degenerate_alphabet_and_complement(self):
        source = '''
        DNA probe = "ATGNNNCRYW";
        PRINT COMPLEMENT(probe);
        '''
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        SemanticAnalyzer().analyze(program)
        tac = TACGenerator().generate(program)
        rt = Runtime(capture_output=True)
        lines = rt.execute(tac)
        # A->T, T->A, G->C, N->N, C->G, R(A/G)->Y(T/C), Y(C/T)->R(G/A), W(A/T)->W(T/A)
        self.assertEqual(lines, ["TACNNNGYRW"])

    def test_fasta_io_load_and_save(self):
        source = '''
        DNA seq = LOAD_FASTA("examples/sample_gene.fasta");
        PRINT LENGTH(seq);
        INTEGER count = SAVE_FASTA("examples/test_out.fasta", seq, "Test Header");
        PRINT count;
        '''
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        SemanticAnalyzer().analyze(program)
        tac = TACGenerator().generate(program)
        rt = Runtime(capture_output=True)
        lines = rt.execute(tac)
        self.assertEqual(lines, ["198", "198"])
        if os.path.exists("examples/test_out.fasta"):
            os.remove("examples/test_out.fasta")

    def test_aho_corasick_fused_motif_search(self):
        source = '''
        DNA seq = "ATGCGTATGCCGATG";
        PRINT COUNT_MOTIF(seq, "ATG");
        PRINT COUNT_MOTIF(seq, "CGT");
        PRINT COUNT_MOTIF(seq, "GGA");
        '''
        tokens = Lexer(source).tokenize()
        program = Parser(tokens).parse_program()
        SemanticAnalyzer().analyze(program)
        tac = TACGenerator().generate(program)
        opt_tac, stats = Optimizer().optimize(tac)
        rt = Runtime(capture_output=True)
        lines = rt.execute(opt_tac)
        self.assertEqual(lines, ["3", "1", "0"])


if __name__ == "__main__":
    unittest.main()
