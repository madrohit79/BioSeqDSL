"""
Lexer for BioSeqDSL.

Implements the lexical-level EBNF rules from the Review 1 document:
IDENTIFIER, STRING_LITERAL, INTEGER_LITERAL, FLOAT_LITERAL, plus
keywords and punctuation. Per the grammar's design note, every quoted
string is lexed as a single STRING_LITERAL — DNA/RNA alphabet
validation happens later, in semantic analysis, once the declared
type is known.
"""

from tokens import Token, TokenType, KEYWORDS
from errors import LexError


class Lexer:
    def __init__(self, source: str):
        self.source = source
        self.pos = 0
        self.line = 1
        self.col = 1
        self.tokens: list[Token] = []

    # ---- low-level cursor helpers ----

    def _peek(self, offset: int = 0) -> str:
        idx = self.pos + offset
        return self.source[idx] if idx < len(self.source) else "\0"

    def _advance(self) -> str:
        ch = self.source[self.pos]
        self.pos += 1
        if ch == "\n":
            self.line += 1
            self.col = 1
        else:
            self.col += 1
        return ch

    def _match(self, expected: str) -> bool:
        if self._peek() != expected:
            return False
        self._advance()
        return True

    def _add_token(self, ttype: TokenType, lexeme: str, line: int, col: int):
        self.tokens.append(Token(ttype, lexeme, line, col))

    # ---- main entry point ----

    def tokenize(self) -> list[Token]:
        while self.pos < len(self.source):
            self._scan_token()
        self.tokens.append(Token(TokenType.EOF, "", self.line, self.col))
        return self.tokens

    def _scan_token(self):
        ch = self._peek()

        # whitespace — skip
        if ch in " \t\r\n":
            self._advance()
            return

        # comments — // to end of line
        if ch == "/" and self._peek(1) == "/":
            while self._peek() not in ("\n", "\0"):
                self._advance()
            return

        start_line, start_col = self.line, self.col

        if ch == '"':
            self._scan_string(start_line, start_col)
            return

        if ch.isdigit():
            self._scan_number(start_line, start_col)
            return

        if ch.isalpha() or ch == "_":
            self._scan_identifier_or_keyword(start_line, start_col)
            return

        # Two-character or single-character operators
        if ch == "=":
            self._advance()
            if self._peek() == "=":
                self._advance()
                self._add_token(TokenType.EQ, "==", start_line, start_col)
            else:
                self._add_token(TokenType.ASSIGN, "=", start_line, start_col)
            return

        if ch == "!":
            self._advance()
            if self._peek() == "=":
                self._advance()
                self._add_token(TokenType.NEQ, "!=", start_line, start_col)
                return
            raise LexError(
                "Unexpected character '!' (did you mean '!='?)",
                start_line,
                start_col,
                hint="Expected '!=' for inequality comparison.",
            )

        if ch == "<":
            self._advance()
            if self._peek() == "=":
                self._advance()
                self._add_token(TokenType.LTE, "<=", start_line, start_col)
            else:
                self._add_token(TokenType.LT, "<", start_line, start_col)
            return

        if ch == ">":
            self._advance()
            if self._peek() == "=":
                self._advance()
                self._add_token(TokenType.GTE, ">=", start_line, start_col)
            else:
                self._add_token(TokenType.GT, ">", start_line, start_col)
            return

        # single-character punctuation
        simple = {
            ";": TokenType.SEMICOLON,
            "(": TokenType.LPAREN,
            ")": TokenType.RPAREN,
            "{": TokenType.LBRACE,
            "}": TokenType.RBRACE,
            ",": TokenType.COMMA,
            "+": TokenType.PLUS,
            "-": TokenType.MINUS,
        }
        if ch in simple:
            self._advance()
            self._add_token(simple[ch], ch, start_line, start_col)
            return

        # anything else is illegal at the lexical level
        self._advance()
        raise LexError(
            f"Unexpected character '{ch}'",
            start_line,
            start_col,
            hint="Check for typos or an unsupported symbol.",
        )

    def _scan_string(self, line: int, col: int):
        self._advance()  # consume opening "
        chars = []
        while self._peek() != '"':
            if self._peek() == "\0":
                raise LexError(
                    "Unterminated string literal",
                    line,
                    col,
                    hint='Did you forget a closing "?',
                )
            chars.append(self._advance())
        self._advance()  # consume closing "
        lexeme = "".join(chars)
        self._add_token(TokenType.STRING_LITERAL, lexeme, line, col)

    def _scan_number(self, line: int, col: int):
        chars = []
        while self._peek().isdigit():
            chars.append(self._advance())

        is_float = False
        if self._peek() == "." and self._peek(1).isdigit():
            is_float = True
            chars.append(self._advance())  # consume '.'
            while self._peek().isdigit():
                chars.append(self._advance())

        lexeme = "".join(chars)
        ttype = TokenType.FLOAT_LITERAL if is_float else TokenType.INTEGER_LITERAL
        self._add_token(ttype, lexeme, line, col)

    def _scan_identifier_or_keyword(self, line: int, col: int):
        chars = []
        while self._peek().isalnum() or self._peek() == "_":
            chars.append(self._advance())
        lexeme = "".join(chars)
        ttype = KEYWORDS.get(lexeme, TokenType.IDENTIFIER)
        self._add_token(ttype, lexeme, line, col)
