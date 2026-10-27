"""
AST node classes for BioSeqDSL.

Every node carries the line/column of the token it was built from,
so the symbol table and later the semantic analyzer can report
exactly where a declaration or error occurred.
"""

from dataclasses import dataclass, field


class ASTNode:
    """Marker base class for all AST nodes."""
    pass


# ---- expressions ----

@dataclass
class Literal(ASTNode):
    value: str          # raw lexeme, e.g. "ATGCGT", "42", "3.14"
    kind: str           # "STRING" | "INTEGER" | "FLOAT"
    line: int
    column: int


@dataclass
class Identifier(ASTNode):
    name: str
    line: int
    column: int


@dataclass
class FunctionCall(ASTNode):
    name: str                      # e.g. "GC_CONTENT"
    arguments: list = field(default_factory=list)   # list[ASTNode]
    line: int = 0
    column: int = 0


@dataclass
class BinaryOp(ASTNode):
    left: ASTNode
    op: str                        # "<", ">", "<=", ">=", "==", "!=", "AND", "OR"
    right: ASTNode
    line: int = 0
    column: int = 0


@dataclass
class UnaryOp(ASTNode):
    op: str                        # "NOT"
    operand: ASTNode
    line: int = 0
    column: int = 0


# ---- statements ----

@dataclass
class Block(ASTNode):
    statements: list = field(default_factory=list)  # list[ASTNode]
    line: int = 0
    column: int = 0


@dataclass
class IfStatement(ASTNode):
    condition: ASTNode
    then_block: Block
    else_block: Block | None = None
    line: int = 0
    column: int = 0


@dataclass
class WhileStatement(ASTNode):
    condition: ASTNode
    body: Block
    line: int = 0
    column: int = 0


@dataclass
class Declaration(ASTNode):
    var_type: str        # "DNA" | "RNA" | "INTEGER" | "FLOAT" | "STRING" | "BOOLEAN"
    name: str
    value: ASTNode        # expression
    line: int
    column: int


@dataclass
class Assignment(ASTNode):
    name: str
    value: ASTNode
    line: int
    column: int


@dataclass
class PrintStatement(ASTNode):
    value: ASTNode
    line: int
    column: int


@dataclass
class Program(ASTNode):
    statements: list = field(default_factory=list)   # list[ASTNode]


def pretty(node: ASTNode, indent: int = 0) -> str:
    """Render an AST as an indented tree — used for the Review 1 demo."""
    pad = "  " * indent

    if isinstance(node, Program):
        lines = [f"{pad}Program"]
        for stmt in node.statements:
            lines.append(pretty(stmt, indent + 1))
        return "\n".join(lines)

    if isinstance(node, Declaration):
        lines = [f"{pad}Declaration({node.var_type} {node.name}) [line {node.line}]"]
        lines.append(pretty(node.value, indent + 1))
        return "\n".join(lines)

    if isinstance(node, Assignment):
        lines = [f"{pad}Assignment({node.name}) [line {node.line}]"]
        lines.append(pretty(node.value, indent + 1))
        return "\n".join(lines)

    if isinstance(node, PrintStatement):
        lines = [f"{pad}Print [line {node.line}]"]
        lines.append(pretty(node.value, indent + 1))
        return "\n".join(lines)

    if isinstance(node, Block):
        lines = [f"{pad}Block [line {node.line}]"]
        for stmt in node.statements:
            lines.append(pretty(stmt, indent + 1))
        return "\n".join(lines)

    if isinstance(node, IfStatement):
        lines = [f"{pad}IfStatement [line {node.line}]"]
        lines.append(f"{pad}  Condition:")
        lines.append(pretty(node.condition, indent + 2))
        lines.append(f"{pad}  Then:")
        lines.append(pretty(node.then_block, indent + 2))
        if node.else_block:
            lines.append(f"{pad}  Else:")
            lines.append(pretty(node.else_block, indent + 2))
        return "\n".join(lines)

    if isinstance(node, WhileStatement):
        lines = [f"{pad}WhileStatement [line {node.line}]"]
        lines.append(f"{pad}  Condition:")
        lines.append(pretty(node.condition, indent + 2))
        lines.append(f"{pad}  Body:")
        lines.append(pretty(node.body, indent + 2))
        return "\n".join(lines)

    if isinstance(node, BinaryOp):
        lines = [f"{pad}BinaryOp({node.op}) [line {node.line}]"]
        lines.append(pretty(node.left, indent + 1))
        lines.append(pretty(node.right, indent + 1))
        return "\n".join(lines)

    if isinstance(node, UnaryOp):
        lines = [f"{pad}UnaryOp({node.op}) [line {node.line}]"]
        lines.append(pretty(node.operand, indent + 1))
        return "\n".join(lines)

    if isinstance(node, FunctionCall):
        lines = [f"{pad}FunctionCall({node.name}) [line {node.line}]"]
        for arg in node.arguments:
            lines.append(pretty(arg, indent + 1))
        return "\n".join(lines)

    if isinstance(node, Literal):
        return f"{pad}Literal({node.kind}, {node.value!r}) [line {node.line}]"

    if isinstance(node, Identifier):
        return f"{pad}Identifier({node.name}) [line {node.line}]"

    return f"{pad}{node!r}"
