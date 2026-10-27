"""
Recursive-descent parser for BioSeqDSL.

Each method below corresponds directly to a production in the
Formal Grammar (EBNF) section of the Review 1 document:

    program          ::= statement*
    statement        ::= declaration | assignment | print_statement
    declaration      ::= type IDENTIFIER '=' expression ';'
    assignment       ::= IDENTIFIER '=' expression ';'
    print_statement  ::= 'PRINT' expression ';'
    expression       ::= literal | IDENTIFIER | function_call
    function_call    ::= FUNCTION_NAME '(' argument_list? ')'
    argument_list    ::= expression (',' expression)*

The grammar is LL(1): parse_statement() always knows which rule to
use from a single token of lookahead (a type keyword -> declaration,
an identifier -> assignment, PRINT -> print_statement), so no
backtracking is needed anywhere in this file.
"""

from tokens import Token, TokenType, TYPE_TOKENS, FUNCTION_TOKENS
from ast_nodes import (
    Program, Declaration, Assignment, PrintStatement,
    FunctionCall, Literal, Identifier, Block,
    IfStatement, WhileStatement, BinaryOp, UnaryOp,
)
from errors import ParseError


class Parser:
    def __init__(self, tokens: list[Token]):
        self.tokens = tokens
        self.pos = 0

    # ---- token stream helpers ----

    def _peek(self) -> Token:
        return self.tokens[self.pos]

    def _advance(self) -> Token:
        tok = self.tokens[self.pos]
        if tok.type != TokenType.EOF:
            self.pos += 1
        return tok

    def _check(self, ttype: TokenType) -> bool:
        return self._peek().type == ttype

    def _expect(self, ttype: TokenType, context: str) -> Token:
        if self._check(ttype):
            return self._advance()
        tok = self._peek()
        raise ParseError(
            f"Expected {ttype.name} {context}, found {tok.type.name} ({tok.lexeme!r})",
            tok.line,
            tok.column,
            hint="Check for a missing token (often ';', '=', '(' or ')').",
        )

    # ---- grammar rules ----

    def parse_program(self) -> Program:
        statements = []
        while not self._check(TokenType.EOF):
            statements.append(self._parse_statement())
        return Program(statements=statements)

    def _parse_statement(self):
        tok = self._peek()

        if tok.type in TYPE_TOKENS:
            return self._parse_declaration()
        if tok.type == TokenType.IDENTIFIER:
            return self._parse_assignment()
        if tok.type == TokenType.PRINT:
            return self._parse_print_statement()
        if tok.type == TokenType.IF:
            return self._parse_if_statement()
        if tok.type == TokenType.WHILE:
            return self._parse_while_statement()
        if tok.type == TokenType.LBRACE:
            return self._parse_block()

        raise ParseError(
            f"Unexpected token {tok.type.name} ({tok.lexeme!r}) at start of statement",
            tok.line,
            tok.column,
            hint="A statement must start with a type keyword, an identifier, PRINT, IF, or WHILE.",
        )

    def _parse_block(self) -> Block:
        lbrace = self._expect(TokenType.LBRACE, "to start code block")
        statements = []
        while not self._check(TokenType.RBRACE) and not self._check(TokenType.EOF):
            statements.append(self._parse_statement())
        self._expect(TokenType.RBRACE, "to close code block")
        return Block(statements=statements, line=lbrace.line, column=lbrace.column)

    def _parse_if_statement(self) -> IfStatement:
        if_tok = self._advance()  # IF
        self._expect(TokenType.LPAREN, "after IF")
        condition = self._parse_expression()
        self._expect(TokenType.RPAREN, "after condition in IF")
        then_block = self._parse_block()
        else_block = None

        if self._check(TokenType.ELSE):
            self._advance()  # ELSE
            if self._check(TokenType.IF):
                # ELSE IF -> nest an IfStatement inside a Block
                nested_if = self._parse_if_statement()
                else_block = Block(statements=[nested_if], line=nested_if.line, column=nested_if.column)
            else:
                else_block = self._parse_block()

        return IfStatement(
            condition=condition,
            then_block=then_block,
            else_block=else_block,
            line=if_tok.line,
            column=if_tok.column,
        )

    def _parse_while_statement(self) -> WhileStatement:
        while_tok = self._advance()  # WHILE
        self._expect(TokenType.LPAREN, "after WHILE")
        condition = self._parse_expression()
        self._expect(TokenType.RPAREN, "after condition in WHILE")
        body = self._parse_block()
        return WhileStatement(
            condition=condition,
            body=body,
            line=while_tok.line,
            column=while_tok.column,
        )

    def _parse_declaration(self) -> Declaration:
        type_tok = self._advance()  # one of TYPE_TOKENS
        name_tok = self._expect(TokenType.IDENTIFIER, "after type in declaration")
        self._expect(TokenType.ASSIGN, "after identifier in declaration")
        value = self._parse_expression()
        self._expect(TokenType.SEMICOLON, "to terminate declaration")
        return Declaration(
            var_type=type_tok.type.name,
            name=name_tok.lexeme,
            value=value,
            line=type_tok.line,
            column=type_tok.column,
        )

    def _parse_assignment(self) -> Assignment:
        name_tok = self._advance()  # IDENTIFIER
        self._expect(TokenType.ASSIGN, "after identifier in assignment")
        value = self._parse_expression()
        self._expect(TokenType.SEMICOLON, "to terminate assignment")
        return Assignment(
            name=name_tok.lexeme,
            value=value,
            line=name_tok.line,
            column=name_tok.column,
        )

    def _parse_print_statement(self) -> PrintStatement:
        print_tok = self._advance()  # PRINT
        value = self._parse_expression()
        self._expect(TokenType.SEMICOLON, "to terminate print statement")
        return PrintStatement(value=value, line=print_tok.line, column=print_tok.column)

    # ---- expression precedence parsing ----

    def _parse_expression(self):
        return self._parse_logical_or()

    def _parse_logical_or(self):
        expr = self._parse_logical_and()
        while self._check(TokenType.OR):
            op_tok = self._advance()
            right = self._parse_logical_and()
            expr = BinaryOp(left=expr, op="OR", right=right, line=op_tok.line, column=op_tok.column)
        return expr

    def _parse_logical_and(self):
        expr = self._parse_comparison()
        while self._check(TokenType.AND):
            op_tok = self._advance()
            right = self._parse_comparison()
            expr = BinaryOp(left=expr, op="AND", right=right, line=op_tok.line, column=op_tok.column)
        return expr

    def _parse_comparison(self):
        expr = self._parse_additive()
        comp_ops = {
            TokenType.EQ: "==",
            TokenType.NEQ: "!=",
            TokenType.LT: "<",
            TokenType.GT: ">",
            TokenType.LTE: "<=",
            TokenType.GTE: ">=",
        }
        while self._peek().type in comp_ops:
            op_tok = self._advance()
            op_symbol = comp_ops[op_tok.type]
            right = self._parse_additive()
            expr = BinaryOp(left=expr, op=op_symbol, right=right, line=op_tok.line, column=op_tok.column)
        return expr

    def _parse_additive(self):
        expr = self._parse_unary()
        while self._peek().type in (TokenType.PLUS, TokenType.MINUS):
            op_tok = self._advance()
            op_symbol = "+" if op_tok.type == TokenType.PLUS else "-"
            right = self._parse_unary()
            expr = BinaryOp(left=expr, op=op_symbol, right=right, line=op_tok.line, column=op_tok.column)
        return expr

    def _parse_unary(self):
        if self._check(TokenType.NOT):
            op_tok = self._advance()
            operand = self._parse_unary()
            return UnaryOp(op="NOT", operand=operand, line=op_tok.line, column=op_tok.column)
        return self._parse_primary()

    def _parse_primary(self):
        tok = self._peek()

        if tok.type == TokenType.LPAREN:
            self._advance()
            expr = self._parse_expression()
            self._expect(TokenType.RPAREN, "after grouped expression")
            return expr

        if tok.type == TokenType.TRUE:
            self._advance()
            return Literal(value="TRUE", kind="BOOLEAN", line=tok.line, column=tok.column)

        if tok.type == TokenType.FALSE:
            self._advance()
            return Literal(value="FALSE", kind="BOOLEAN", line=tok.line, column=tok.column)

        if tok.type == TokenType.STRING_LITERAL:
            self._advance()
            return Literal(value=tok.lexeme, kind="STRING", line=tok.line, column=tok.column)

        if tok.type == TokenType.INTEGER_LITERAL:
            self._advance()
            return Literal(value=tok.lexeme, kind="INTEGER", line=tok.line, column=tok.column)

        if tok.type == TokenType.FLOAT_LITERAL:
            self._advance()
            return Literal(value=tok.lexeme, kind="FLOAT", line=tok.line, column=tok.column)

        if tok.type in FUNCTION_TOKENS:
            return self._parse_function_call()

        if tok.type == TokenType.IDENTIFIER:
            self._advance()
            return Identifier(name=tok.lexeme, line=tok.line, column=tok.column)

        raise ParseError(
            f"Expected an expression, found {tok.type.name} ({tok.lexeme!r})",
            tok.line,
            tok.column,
            hint="Expected a literal, identifier, condition, or function call here.",
        )

    def _parse_function_call(self) -> FunctionCall:
        name_tok = self._advance()  # one of FUNCTION_TOKENS
        self._expect(TokenType.LPAREN, "after function name")

        arguments = []
        if not self._check(TokenType.RPAREN):
            arguments.append(self._parse_expression())
            while self._check(TokenType.COMMA):
                self._advance()
                arguments.append(self._parse_expression())

        self._expect(TokenType.RPAREN, "to close function call")
        return FunctionCall(
            name=name_tok.lexeme,
            arguments=arguments,
            line=name_tok.line,
            column=name_tok.column,
        )
