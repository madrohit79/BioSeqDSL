"""
Symbol table for BioSeqDSL.

Walks the AST produced by the parser and records, for every declared
identifier: its declared type and the line it was declared on. Also
catches the two structural semantic errors that don't require full
type checking: using an identifier that was never declared, and
declaring the same identifier twice.

Full domain-aware type checking (e.g. TRANSLATE requires RNA) is the
next phase (semantic_analyzer.py) and builds on top of this table.
"""

from dataclasses import dataclass

from ast_nodes import (
    Program, Declaration, Assignment, PrintStatement,
    FunctionCall, Literal, Identifier, Block,
    IfStatement, WhileStatement, BinaryOp, UnaryOp, ASTNode,
)
from errors import UndeclaredIdentifierError, RedeclarationError


@dataclass
class SymbolInfo:
    var_type: str
    declared_line: int


class SymbolTable:
    def __init__(self):
        self._symbols: dict[str, SymbolInfo] = {}

    def declare(self, name: str, var_type: str, line: int, column: int):
        if name in self._symbols:
            existing = self._symbols[name]
            raise RedeclarationError(
                f"'{name}' is already declared as {existing.var_type} on line {existing.declared_line}",
                line,
                column,
                hint=f"Choose a different name, or remove the duplicate declaration of '{name}'.",
            )
        self._symbols[name] = SymbolInfo(var_type=var_type, declared_line=line)

    def lookup(self, name: str, line: int, column: int) -> SymbolInfo:
        if name not in self._symbols:
            raise UndeclaredIdentifierError(
                f"'{name}' is used but was never declared",
                line,
                column,
                hint=f"Declare it first, e.g. DNA {name} = \"...\";",
            )
        return self._symbols[name]

    def __repr__(self) -> str:
        if not self._symbols:
            return "SymbolTable (empty)"
        lines = ["SymbolTable:"]
        for name, info in self._symbols.items():
            lines.append(f"  {name:<12} type={info.var_type:<8} declared at line {info.declared_line}")
        return "\n".join(lines)


def _check_expression_identifiers(node: ASTNode, table: SymbolTable):
    """Recursively verify every identifier referenced inside an expression is declared."""
    if isinstance(node, Identifier):
        table.lookup(node.name, node.line, node.column)
    elif isinstance(node, FunctionCall):
        for arg in node.arguments:
            _check_expression_identifiers(arg, table)
    elif isinstance(node, BinaryOp):
        _check_expression_identifiers(node.left, table)
        _check_expression_identifiers(node.right, table)
    elif isinstance(node, UnaryOp):
        _check_expression_identifiers(node.operand, table)
    elif isinstance(node, Literal):
        pass  # nothing to check
    else:
        raise TypeError(f"Unexpected expression node: {node!r}")


def _walk_statements(stmts: list[ASTNode], table: SymbolTable):
    for stmt in stmts:
        if isinstance(stmt, Declaration):
            _check_expression_identifiers(stmt.value, table)
            table.declare(stmt.name, stmt.var_type, stmt.line, stmt.column)

        elif isinstance(stmt, Assignment):
            table.lookup(stmt.name, stmt.line, stmt.column)
            _check_expression_identifiers(stmt.value, table)

        elif isinstance(stmt, PrintStatement):
            _check_expression_identifiers(stmt.value, table)

        elif isinstance(stmt, Block):
            _walk_statements(stmt.statements, table)

        elif isinstance(stmt, IfStatement):
            _check_expression_identifiers(stmt.condition, table)
            _walk_statements(stmt.then_block.statements, table)
            if stmt.else_block:
                _walk_statements(stmt.else_block.statements, table)

        elif isinstance(stmt, WhileStatement):
            _check_expression_identifiers(stmt.condition, table)
            _walk_statements(stmt.body.statements, table)

        else:
            raise TypeError(f"Unexpected statement node: {stmt!r}")


def build_symbol_table(program: Program) -> SymbolTable:
    """Walk the AST, populate and return the symbol table.

    Raises RedeclarationError / UndeclaredIdentifierError (both
    CompilerError subclasses) on the first problem found.
    """
    table = SymbolTable()
    _walk_statements(program.statements, table)
    return table
