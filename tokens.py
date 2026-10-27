"""
Token definitions for BioSeqDSL.

A Token carries its line/column so that every later stage (parser,
semantic analyzer) can report precise diagnostics without re-scanning
the source.
"""

from enum import Enum, auto
from dataclasses import dataclass


class TokenType(Enum):
    # Type keywords
    DNA = auto()
    RNA = auto()
    INTEGER = auto()
    FLOAT = auto()
    STRING = auto()
    BOOLEAN = auto()

    # Statement keywords
    PRINT = auto()
    IF = auto()
    ELSE = auto()
    WHILE = auto()

    # Boolean literals
    TRUE = auto()
    FALSE = auto()

    # Logical operators
    AND = auto()
    OR = auto()
    NOT = auto()

    # Built-in function names
    LENGTH = auto()
    GC_CONTENT = auto()
    REVERSE = auto()
    COMPLEMENT = auto()
    REVERSE_COMPLEMENT = auto()
    TRANSCRIBE = auto()
    TRANSLATE = auto()
    FIND_MOTIF = auto()
    COUNT_MOTIF = auto()
    LOAD_FASTA = auto()
    SAVE_FASTA = auto()

    # Literals
    IDENTIFIER = auto()
    STRING_LITERAL = auto()
    INTEGER_LITERAL = auto()
    FLOAT_LITERAL = auto()

    # Punctuation / operators
    ASSIGN = auto()       # =
    SEMICOLON = auto()    # ;
    LPAREN = auto()       # (
    RPAREN = auto()       # )
    LBRACE = auto()       # {
    RBRACE = auto()       # }
    COMMA = auto()        # ,

    # Comparison operators
    EQ = auto()           # ==
    NEQ = auto()          # !=
    LTE = auto()          # <=
    GTE = auto()          # >=
    LT = auto()           # <
    GT = auto()           # >

    # Arithmetic operators
    PLUS = auto()         # +
    MINUS = auto()        # -

    EOF = auto()


# Reserved words -> token type. Checked after an identifier is scanned.
KEYWORDS = {
    "DNA": TokenType.DNA,
    "RNA": TokenType.RNA,
    "INTEGER": TokenType.INTEGER,
    "FLOAT": TokenType.FLOAT,
    "STRING": TokenType.STRING,
    "BOOLEAN": TokenType.BOOLEAN,
    "PRINT": TokenType.PRINT,
    "IF": TokenType.IF,
    "ELSE": TokenType.ELSE,
    "WHILE": TokenType.WHILE,
    "TRUE": TokenType.TRUE,
    "FALSE": TokenType.FALSE,
    "AND": TokenType.AND,
    "OR": TokenType.OR,
    "NOT": TokenType.NOT,
    "LENGTH": TokenType.LENGTH,
    "GC_CONTENT": TokenType.GC_CONTENT,
    "REVERSE": TokenType.REVERSE,
    "COMPLEMENT": TokenType.COMPLEMENT,
    "REVERSE_COMPLEMENT": TokenType.REVERSE_COMPLEMENT,
    "TRANSCRIBE": TokenType.TRANSCRIBE,
    "TRANSLATE": TokenType.TRANSLATE,
    "FIND_MOTIF": TokenType.FIND_MOTIF,
    "COUNT_MOTIF": TokenType.COUNT_MOTIF,
    "LOAD_FASTA": TokenType.LOAD_FASTA,
    "SAVE_FASTA": TokenType.SAVE_FASTA,
}

# The subset of keywords that are valid builtin function names —
# used by the parser to recognize function_call vs plain identifier.
FUNCTION_TOKENS = {
    TokenType.LENGTH,
    TokenType.GC_CONTENT,
    TokenType.REVERSE,
    TokenType.COMPLEMENT,
    TokenType.REVERSE_COMPLEMENT,
    TokenType.TRANSCRIBE,
    TokenType.TRANSLATE,
    TokenType.FIND_MOTIF,
    TokenType.COUNT_MOTIF,
    TokenType.LOAD_FASTA,
    TokenType.SAVE_FASTA,
}

# The subset of keywords that are valid type names — used by the
# parser to recognize the start of a declaration.
TYPE_TOKENS = {
    TokenType.DNA,
    TokenType.RNA,
    TokenType.INTEGER,
    TokenType.FLOAT,
    TokenType.STRING,
    TokenType.BOOLEAN,
}


@dataclass
class Token:
    type: TokenType
    lexeme: str
    line: int
    column: int

    def __repr__(self) -> str:
        return f"Token({self.type.name}, {self.lexeme!r}, {self.line}:{self.column})"
