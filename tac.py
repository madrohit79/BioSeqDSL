"""
Intermediate Code Generation (Three-Address Code / TAC) for BioSeqDSL.

Lowers the checked Abstract Syntax Tree into linear Three-Address Code
instructions (Module 4 of the compiler syllabus).

Supported TAC instructions:
  - ASSIGN: dest = src
  - CALL:   dest = CALL func_name, (args...)
  - PRINT:  PRINT arg
  - FUSED_COUNT_MOTIF: [dest1, dest2, ...] = FUSED_COUNT_MOTIF seq, [pat1, pat2, ...]
"""

from dataclasses import dataclass, field
from ast_nodes import (
    Program, Declaration, Assignment, PrintStatement,
    FunctionCall, Literal, Identifier, Block,
    IfStatement, WhileStatement, BinaryOp, UnaryOp, ASTNode,
)


@dataclass
class TACInstruction:
    op: str                      # 'ASSIGN', 'CALL', 'PRINT', 'LABEL', 'GOTO', 'IF_FALSE', 'BINARY_OP', 'UNARY_OP', 'FUSED_COUNT_MOTIF', 'FUSED_FIND_MOTIF'
    dest: str | list[str] | None # destination variable/temp, or list for fused ops
    args: list[str] = field(default_factory=list) # operands / arguments
    comment: str | None = None
    line: int = 0

    def __repr__(self) -> str:
        if self.op == "ASSIGN":
            return f"{self.dest} = {self.args[0]}"
        if self.op == "CALL":
            func_name = self.args[0]
            call_args = ", ".join(self.args[1:])
            return f"{self.dest} = CALL {func_name}({call_args})"
        if self.op == "PRINT":
            return f"PRINT {self.args[0]}"
        if self.op == "LABEL":
            return f"LABEL {self.args[0]}:"
        if self.op == "GOTO":
            return f"GOTO {self.args[0]}"
        if self.op == "IF_FALSE":
            return f"IF_FALSE {self.args[0]} GOTO {self.args[1]}"
        if self.op == "BINARY_OP":
            return f"{self.dest} = {self.args[1]} {self.args[0]} {self.args[2]}"
        if self.op == "UNARY_OP":
            return f"{self.dest} = {self.args[0]} {self.args[1]}"
        if self.op == "FUSED_COUNT_MOTIF":
            dests = ", ".join(self.dest) if isinstance(self.dest, list) else str(self.dest)
            seq = self.args[0]
            patterns = ", ".join(repr(p) for p in self.args[1:])
            return f"[{dests}] = FUSED_COUNT_MOTIF {seq}, [{patterns}]"
        if self.op == "FUSED_FIND_MOTIF":
            dests = ", ".join(self.dest) if isinstance(self.dest, list) else str(self.dest)
            seq = self.args[0]
            patterns = ", ".join(repr(p) for p in self.args[1:])
            return f"[{dests}] = FUSED_FIND_MOTIF {seq}, [{patterns}]"
        return f"{self.op} {self.dest} {self.args}"


class TACProgram:
    def __init__(self):
        self.instructions: list[TACInstruction] = []

    def emit(self, instruction: TACInstruction):
        self.instructions.append(instruction)

    def render(self) -> str:
        lines = []
        for i, inst in enumerate(self.instructions, 1):
            comment_str = f"    ; {inst.comment}" if inst.comment else ""
            lines.append(f"{i:03d}:  {repr(inst)}{comment_str}")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.render()


class TACGenerator:
    def __init__(self):
        self.temp_count = 0
        self.label_count = 0
        self.tac = TACProgram()

    def _new_temp(self) -> str:
        self.temp_count += 1
        return f"t{self.temp_count}"

    def _new_label(self) -> str:
        self.label_count += 1
        return f"L{self.label_count}"

    def generate(self, program: Program) -> TACProgram:
        for stmt in program.statements:
            self._lower_statement(stmt)
        return self.tac

    def _lower_statement(self, stmt: ASTNode):
        if isinstance(stmt, Declaration):
            rhs_val = self._lower_expression(stmt.value)
            self.tac.emit(TACInstruction(
                op="ASSIGN",
                dest=stmt.name,
                args=[rhs_val],
                comment=f"declare {stmt.var_type} {stmt.name}",
                line=stmt.line,
            ))
        elif isinstance(stmt, Assignment):
            rhs_val = self._lower_expression(stmt.value)
            self.tac.emit(TACInstruction(
                op="ASSIGN",
                dest=stmt.name,
                args=[rhs_val],
                comment=f"assign to {stmt.name}",
                line=stmt.line,
            ))
        elif isinstance(stmt, PrintStatement):
            val = self._lower_expression(stmt.value)
            self.tac.emit(TACInstruction(
                op="PRINT",
                dest=None,
                args=[val],
                line=stmt.line,
            ))
        elif isinstance(stmt, Block):
            for s in stmt.statements:
                self._lower_statement(s)
        elif isinstance(stmt, IfStatement):
            cond_temp = self._lower_expression(stmt.condition)
            else_label = self._new_label()
            end_label = self._new_label()
            self.tac.emit(TACInstruction(
                op="IF_FALSE",
                dest=None,
                args=[cond_temp, else_label],
                line=stmt.line,
            ))
            self._lower_statement(stmt.then_block)
            if stmt.else_block:
                self.tac.emit(TACInstruction(
                    op="GOTO",
                    dest=None,
                    args=[end_label],
                    line=stmt.line,
                ))
            self.tac.emit(TACInstruction(
                op="LABEL",
                dest=None,
                args=[else_label],
                line=stmt.line,
            ))
            if stmt.else_block:
                self._lower_statement(stmt.else_block)
                self.tac.emit(TACInstruction(
                    op="LABEL",
                    dest=None,
                    args=[end_label],
                    line=stmt.line,
                ))
        elif isinstance(stmt, WhileStatement):
            start_label = self._new_label()
            end_label = self._new_label()
            self.tac.emit(TACInstruction(
                op="LABEL",
                dest=None,
                args=[start_label],
                line=stmt.line,
            ))
            cond_temp = self._lower_expression(stmt.condition)
            self.tac.emit(TACInstruction(
                op="IF_FALSE",
                dest=None,
                args=[cond_temp, end_label],
                line=stmt.line,
            ))
            self._lower_statement(stmt.body)
            self.tac.emit(TACInstruction(
                op="GOTO",
                dest=None,
                args=[start_label],
                line=stmt.line,
            ))
            self.tac.emit(TACInstruction(
                op="LABEL",
                dest=None,
                args=[end_label],
                line=stmt.line,
            ))
        else:
            raise TypeError(f"Cannot lower statement: {stmt!r}")

    def _lower_expression(self, expr: ASTNode) -> str:
        if isinstance(expr, Literal):
            if expr.kind == "STRING":
                return f'"{expr.value}"'
            return expr.value

        if isinstance(expr, Identifier):
            return expr.name

        if isinstance(expr, FunctionCall):
            lowered_args = [self._lower_expression(arg) for arg in expr.arguments]
            t = self._new_temp()
            self.tac.emit(TACInstruction(
                op="CALL",
                dest=t,
                args=[expr.name] + lowered_args,
                line=expr.line,
            ))
            return t

        if isinstance(expr, BinaryOp):
            left_val = self._lower_expression(expr.left)
            right_val = self._lower_expression(expr.right)
            t = self._new_temp()
            self.tac.emit(TACInstruction(
                op="BINARY_OP",
                dest=t,
                args=[expr.op, left_val, right_val],
                line=expr.line,
            ))
            return t

        if isinstance(expr, UnaryOp):
            val = self._lower_expression(expr.operand)
            t = self._new_temp()
            self.tac.emit(TACInstruction(
                op="UNARY_OP",
                dest=t,
                args=[expr.op, val],
                line=expr.line,
            ))
            return t

        raise TypeError(f"Cannot lower expression: {expr!r}")
