"""
Shared diagnostic machinery for BioSeqDSL.

Every stage (lexer, parser, semantic analyzer) raises a CompilerError
subclass instead of crashing. main.py catches these and prints them
in the line/column + category + snippet + hint format described in
the Review 1 document's "Semantic Error Diagnostics" section.
"""


class CompilerError(Exception):
    """Base class for all BioSeqDSL compiler diagnostics."""

    category = "Error"

    def __init__(self, message: str, line: int, column: int, hint: str | None = None, span_length: int = 1):
        super().__init__(message)
        self.message = message
        self.line = line
        self.column = column
        self.hint = hint
        self.span_length = max(1, span_length)

    def render(self, source_lines: list[str]) -> str:
        """Format this error the way it's shown in the project doc."""
        out = [f"Error: {self.category} (line {self.line}, column {self.column})"]

        if 1 <= self.line <= len(source_lines):
            src_line = source_lines[self.line - 1]
            out.append(f"  {src_line}")
            caret_pos = max(self.column - 1, 0)
            out.append("  " + " " * caret_pos + "^" * self.span_length)

        out.append(f"  {self.message}")
        if self.hint:
            out.append(f"  Hint: {self.hint}")
        return "\n".join(out)


class LexError(CompilerError):
    category = "Lexical Error"


class ParseError(CompilerError):
    category = "Syntax Error"


class UndeclaredIdentifierError(CompilerError):
    category = "Undeclared Identifier"


class RedeclarationError(CompilerError):
    category = "Redeclaration"


class TypeMismatchError(CompilerError):
    category = "Type Mismatch"


class InvalidArgumentTypeError(CompilerError):
    category = "Invalid Argument Type"


class InvalidAlphabetError(CompilerError):
    category = "Invalid Sequence Alphabet"

