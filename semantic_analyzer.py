"""
Semantic Analyzer for BioSeqDSL.

Performs domain-aware biological semantic checks and type analysis:
1. Validates sequence alphabets for DNA ({A, T, G, C}) and RNA ({A, U, G, C}).
2. Verifies operation signatures (e.g. TRANSLATE requires RNA, TRANSCRIBE requires DNA).
3. Enforces type consistency on declarations and assignments.
4. Generates rich, actionable compiler diagnostics with source spans and hints.
"""

from ast_nodes import (
    Program, Declaration, Assignment, PrintStatement,
    FunctionCall, Literal, Identifier, Block,
    IfStatement, WhileStatement, BinaryOp, UnaryOp, ASTNode,
)
from symbol_table import SymbolTable
from errors import (
    TypeMismatchError, InvalidArgumentTypeError, InvalidAlphabetError,
)

DNA_STANDARD_ALPHABET = set("ATGCatgc")
DNA_IUPAC_ALPHABET = set("ATGCRYSWKMBDHVNatgcryswkmbdhvn")
RNA_STANDARD_ALPHABET = set("AUGCaugc")
RNA_IUPAC_ALPHABET = set("AUGCRYSWKMBDHVNaugcryswkmbdhvn")

DNA_ALPHABET = DNA_IUPAC_ALPHABET
RNA_ALPHABET = RNA_IUPAC_ALPHABET


class SemanticAnalyzer:
    def __init__(self, table: SymbolTable | None = None):
        self.table = table if table is not None else SymbolTable()
        self.node_types: dict[int, str] = {}  # id(node) -> inferred type

    def analyze(self, program: Program) -> SymbolTable:
        """Walk the AST, check types and sequence semantics, return populated symbol table."""
        for stmt in program.statements:
            self._check_statement(stmt)
        return self.table

    def _check_statement(self, stmt: ASTNode):
        if isinstance(stmt, Declaration):
            self._check_declaration(stmt)
        elif isinstance(stmt, Assignment):
            self._check_assignment(stmt)
        elif isinstance(stmt, PrintStatement):
            self._check_print(stmt)
        elif isinstance(stmt, Block):
            self._check_block(stmt)
        elif isinstance(stmt, IfStatement):
            self._check_if(stmt)
        elif isinstance(stmt, WhileStatement):
            self._check_while(stmt)
        else:
            raise TypeError(f"Unknown statement node: {stmt!r}")

    def _check_block(self, block: Block):
        for stmt in block.statements:
            self._check_statement(stmt)

    def _check_if(self, stmt: IfStatement):
        cond_type = self._check_expression(stmt.condition)
        if cond_type != "BOOLEAN":
            raise TypeMismatchError(
                f"Expected: BOOLEAN Found: {cond_type}",
                stmt.condition.line,
                stmt.condition.column,
                hint="Condition in IF statement must evaluate to BOOLEAN.",
                span_length=2,
            )
        self._check_block(stmt.then_block)
        if stmt.else_block:
            self._check_block(stmt.else_block)

    def _check_while(self, stmt: WhileStatement):
        cond_type = self._check_expression(stmt.condition)
        if cond_type != "BOOLEAN":
            raise TypeMismatchError(
                f"Expected: BOOLEAN Found: {cond_type}",
                stmt.condition.line,
                stmt.condition.column,
                hint="Condition in WHILE loop must evaluate to BOOLEAN.",
                span_length=5,
            )
        self._check_block(stmt.body)

    def _check_declaration(self, decl: Declaration):
        val_type = self._check_expression(decl.value)

        # Literal alphabet validation
        if isinstance(decl.value, Literal) and decl.value.kind == "STRING":
            literal_val = decl.value.value
            if decl.var_type == "DNA":
                for ch in literal_val:
                    if ch not in DNA_ALPHABET:
                        raise InvalidAlphabetError(
                            f"Invalid character '{ch}' in DNA sequence literal. DNA only allows A, T, G, C.",
                            decl.value.line,
                            decl.value.column,
                            hint="Check for typos or RNA bases ('U'). DNA nucleotides must be A, T, G, or C.",
                            span_length=len(literal_val) + 2,
                        )
                val_type = "DNA"
            elif decl.var_type == "RNA":
                for ch in literal_val:
                    if ch not in RNA_ALPHABET:
                        raise InvalidAlphabetError(
                            f"Invalid character '{ch}' in RNA sequence literal. RNA uses Uracil (U) instead of Thymine (T).",
                            decl.value.line,
                            decl.value.column,
                            hint="Did you mean to transcribe a DNA sequence, or replace 'T' with 'U'?",
                            span_length=len(literal_val) + 2,
                        )
                val_type = "RNA"

        # Type compatibility check
        if val_type != decl.var_type:
            span_len = len(decl.value.name) if isinstance(decl.value, Identifier) else 3
            hint = f"did you mean to declare '{decl.name}' as {val_type}, or convert first?"
            raise TypeMismatchError(
                f"Expected: {decl.var_type} Found: {val_type}",
                decl.value.line,
                decl.value.column,
                hint=hint,
                span_length=span_len,
            )

        self.table.declare(decl.name, decl.var_type, decl.line, decl.column)
        self.node_types[id(decl)] = decl.var_type

    def _check_assignment(self, assign: Assignment):
        info = self.table.lookup(assign.name, assign.line, assign.column)
        val_type = self._check_expression(assign.value)

        if val_type != info.var_type:
            span_len = len(assign.value.name) if isinstance(assign.value, Identifier) else 3
            raise TypeMismatchError(
                f"Expected: {info.var_type} Found: {val_type}",
                assign.value.line,
                assign.value.column,
                hint=f"Variable '{assign.name}' was declared as {info.var_type} on line {info.declared_line}.",
                span_length=span_len,
            )
        self.node_types[id(assign)] = info.var_type

    def _check_print(self, stmt: PrintStatement):
        val_type = self._check_expression(stmt.value)
        self.node_types[id(stmt)] = val_type

    def _check_expression(self, expr: ASTNode) -> str:
        """Infer and return the type of an expression AST node."""
        if isinstance(expr, Literal):
            # Base literal types
            if expr.kind == "INTEGER":
                res = "INTEGER"
            elif expr.kind == "FLOAT":
                res = "FLOAT"
            elif expr.kind == "STRING":
                res = "STRING"
            else:
                res = expr.kind
            self.node_types[id(expr)] = res
            return res

        if isinstance(expr, Identifier):
            info = self.table.lookup(expr.name, expr.line, expr.column)
            self.node_types[id(expr)] = info.var_type
            return info.var_type

        if isinstance(expr, FunctionCall):
            res = self._check_function_call(expr)
            self.node_types[id(expr)] = res
            return res

        if isinstance(expr, BinaryOp):
            left_t = self._check_expression(expr.left)
            right_t = self._check_expression(expr.right)
            if expr.op in ("<", ">", "<=", ">="):
                if left_t in ("INTEGER", "FLOAT") and right_t in ("INTEGER", "FLOAT"):
                    res = "BOOLEAN"
                else:
                    raise TypeMismatchError(
                        f"Comparison '{expr.op}' requires numeric types, found {left_t} and {right_t}",
                        expr.line,
                        expr.column,
                        hint="Relational comparisons operate on INTEGER or FLOAT.",
                        span_length=len(expr.op),
                    )
            elif expr.op in ("==", "!="):
                res = "BOOLEAN"
            elif expr.op in ("+", "-"):
                if left_t in ("INTEGER", "FLOAT") and right_t in ("INTEGER", "FLOAT"):
                    res = "INTEGER" if left_t == "INTEGER" and right_t == "INTEGER" else "FLOAT"
                elif expr.op == "+" and (left_t in ("STRING", "DNA", "RNA") or right_t in ("STRING", "DNA", "RNA")):
                    res = left_t
                else:
                    raise TypeMismatchError(
                        f"Operator '{expr.op}' requires numeric or sequence types, found {left_t} and {right_t}",
                        expr.line,
                        expr.column,
                        hint="Numeric arithmetic works on INTEGER/FLOAT; '+' concatenates strings or sequences.",
                        span_length=len(expr.op),
                    )
            elif expr.op in ("AND", "OR"):
                if left_t == "BOOLEAN" and right_t == "BOOLEAN":
                    res = "BOOLEAN"
                else:
                    raise TypeMismatchError(
                        f"Logical operator '{expr.op}' requires BOOLEAN operands, found {left_t} and {right_t}",
                        expr.line,
                        expr.column,
                        hint="Logical AND / OR can only combine boolean conditions.",
                        span_length=len(expr.op),
                    )
            else:
                res = "BOOLEAN"
            self.node_types[id(expr)] = res
            return res

        if isinstance(expr, UnaryOp):
            operand_t = self._check_expression(expr.operand)
            if expr.op == "NOT":
                if operand_t != "BOOLEAN":
                    raise TypeMismatchError(
                        f"NOT expects BOOLEAN, received {operand_t}",
                        expr.line,
                        expr.column,
                        hint="NOT operator operates only on boolean conditions.",
                        span_length=3,
                    )
                res = "BOOLEAN"
            else:
                res = operand_t
            self.node_types[id(expr)] = res
            return res

        raise TypeError(f"Unexpected expression node: {expr!r}")

    def _check_function_call(self, call: FunctionCall) -> str:
        name = call.name
        args = call.arguments
        arg_types = [self._check_expression(a) for a in args]

        # 1. LENGTH(seq) -> INTEGER
        if name == "LENGTH":
            self._expect_arity(call, 1, len(args))
            seq_type = arg_types[0]
            if seq_type not in ("DNA", "RNA", "STRING"):
                raise InvalidArgumentTypeError(
                    f"LENGTH expects DNA or RNA, received {seq_type}",
                    call.line,
                    call.column,
                    hint="LENGTH can only be computed on sequence or string types.",
                    span_length=len(name),
                )
            return "INTEGER"

        # 2. GC_CONTENT(seq) -> FLOAT
        if name == "GC_CONTENT":
            self._expect_arity(call, 1, len(args))
            seq_type = arg_types[0]
            if seq_type not in ("DNA", "RNA"):
                raise InvalidArgumentTypeError(
                    f"GC_CONTENT expects DNA or RNA, received {seq_type}",
                    call.line,
                    call.column,
                    hint="GC_CONTENT can only be calculated for nucleotide sequences (DNA or RNA).",
                    span_length=len(name),
                )
            return "FLOAT"

        # 3. REVERSE(seq) -> same type
        if name == "REVERSE":
            self._expect_arity(call, 1, len(args))
            seq_type = arg_types[0]
            if seq_type not in ("DNA", "RNA", "STRING"):
                raise InvalidArgumentTypeError(
                    f"REVERSE expects DNA or RNA, received {seq_type}",
                    call.line,
                    call.column,
                    hint="REVERSE operates on DNA or RNA sequences.",
                    span_length=len(name),
                )
            return seq_type

        # 4. COMPLEMENT(seq) -> same type
        if name == "COMPLEMENT":
            self._expect_arity(call, 1, len(args))
            seq_type = arg_types[0]
            if seq_type not in ("DNA", "RNA"):
                raise InvalidArgumentTypeError(
                    f"COMPLEMENT expects DNA or RNA, received {seq_type}",
                    call.line,
                    call.column,
                    hint="COMPLEMENT operates only on nucleotide sequences (DNA or RNA).",
                    span_length=len(name),
                )
            return seq_type

        # 5. REVERSE_COMPLEMENT(seq) -> same type
        if name == "REVERSE_COMPLEMENT":
            self._expect_arity(call, 1, len(args))
            seq_type = arg_types[0]
            if seq_type not in ("DNA", "RNA"):
                raise InvalidArgumentTypeError(
                    f"REVERSE_COMPLEMENT expects DNA or RNA, received {seq_type}",
                    call.line,
                    call.column,
                    hint="REVERSE_COMPLEMENT operates on nucleotide sequences (DNA or RNA).",
                    span_length=len(name),
                )
            return seq_type

        # 6. TRANSCRIBE(seq) -> RNA
        if name == "TRANSCRIBE":
            self._expect_arity(call, 1, len(args))
            seq_type = arg_types[0]
            if seq_type != "DNA":
                raise InvalidArgumentTypeError(
                    f"TRANSCRIBE expects DNA, received {seq_type}",
                    call.line,
                    call.column,
                    hint="TRANSCRIBE synthesizes RNA from a DNA template. Did you pass an RNA variable by mistake?",
                    span_length=len(name),
                )
            return "RNA"

        # 7. TRANSLATE(seq) -> STRING (protein)
        if name == "TRANSLATE":
            self._expect_arity(call, 1, len(args))
            seq_type = arg_types[0]
            if seq_type != "RNA":
                raise InvalidArgumentTypeError(
                    f"TRANSLATE expects RNA, received {seq_type}",
                    call.line,
                    call.column,
                    hint="TRANSLATE operates on RNA. Did you forget TRANSCRIBE(seq) first?",
                    span_length=len(name),
                )
            return "STRING"

        # 8. FIND_MOTIF(seq, pattern) -> LIST[INTEGER]
        if name == "FIND_MOTIF":
            self._expect_arity(call, 2, len(args))
            seq_type, pat_type = arg_types[0], arg_types[1]
            if seq_type not in ("DNA", "RNA"):
                raise InvalidArgumentTypeError(
                    f"FIND_MOTIF target sequence must be DNA or RNA, received {seq_type}",
                    call.line,
                    call.column,
                    hint="The first argument to FIND_MOTIF must be a sequence variable.",
                    span_length=len(name),
                )
            if pat_type not in ("STRING", "DNA", "RNA"):
                raise InvalidArgumentTypeError(
                    f"FIND_MOTIF motif pattern must be a string or sequence, received {pat_type}",
                    call.line,
                    call.column,
                    hint="Pass the motif pattern as a string literal, e.g. \"ATG\".",
                    span_length=len(name),
                )
            return "LIST[INTEGER]"

        # 9. COUNT_MOTIF(seq, pattern) -> INTEGER
        if name == "COUNT_MOTIF":
            self._expect_arity(call, 2, len(args))
            seq_type, pat_type = arg_types[0], arg_types[1]
            if seq_type not in ("DNA", "RNA"):
                raise InvalidArgumentTypeError(
                    f"COUNT_MOTIF target sequence must be DNA or RNA, received {seq_type}",
                    call.line,
                    call.column,
                    hint="The first argument to COUNT_MOTIF must be a sequence variable.",
                    span_length=len(name),
                )
            if pat_type not in ("STRING", "DNA", "RNA"):
                raise InvalidArgumentTypeError(
                    f"COUNT_MOTIF motif pattern must be a string or sequence, received {pat_type}",
                    call.line,
                    call.column,
                    hint="Pass the motif pattern as a string literal, e.g. \"ATG\".",
                    span_length=len(name),
                )
            return "INTEGER"

        # 10. LOAD_FASTA(filepath) -> DNA
        if name == "LOAD_FASTA":
            self._expect_arity(call, 1, len(args))
            path_type = arg_types[0]
            if path_type != "STRING":
                raise InvalidArgumentTypeError(
                    f"LOAD_FASTA expects a STRING file path, received {path_type}",
                    call.line,
                    call.column,
                    hint="Provide a FASTA file path as a string literal, e.g. LOAD_FASTA(\"genome.fasta\").",
                    span_length=len(name),
                )
            return "DNA"

        # 11. SAVE_FASTA(filepath, seq, header) -> INTEGER
        if name == "SAVE_FASTA":
            self._expect_arity(call, 3, len(args))
            path_type, seq_type, hdr_type = arg_types[0], arg_types[1], arg_types[2]
            if path_type != "STRING":
                raise InvalidArgumentTypeError(
                    f"SAVE_FASTA argument 1 (path) must be STRING, received {path_type}",
                    call.line, call.column, span_length=len(name),
                )
            if seq_type not in ("DNA", "RNA", "STRING"):
                raise InvalidArgumentTypeError(
                    f"SAVE_FASTA argument 2 (sequence) must be DNA, RNA, or STRING, received {seq_type}",
                    call.line, call.column, span_length=len(name),
                )
            if hdr_type != "STRING":
                raise InvalidArgumentTypeError(
                    f"SAVE_FASTA argument 3 (header) must be STRING, received {hdr_type}",
                    call.line, call.column, span_length=len(name),
                )
            return "INTEGER"

        raise InvalidArgumentTypeError(
            f"Unknown function '{name}'",
            call.line,
            call.column,
            hint="Valid functions are LENGTH, GC_CONTENT, REVERSE, COMPLEMENT, REVERSE_COMPLEMENT, TRANSCRIBE, TRANSLATE, FIND_MOTIF, COUNT_MOTIF, LOAD_FASTA, SAVE_FASTA.",
            span_length=len(name),
        )

    def _expect_arity(self, call: FunctionCall, expected: int, found: int):
        if found != expected:
            raise InvalidArgumentTypeError(
                f"{call.name} takes {expected} argument(s), found {found}",
                call.line,
                call.column,
                hint=f"Check the signature of {call.name}(...).",
                span_length=len(call.name),
            )
