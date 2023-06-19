%{
#include <iostream>
#include <string>
#include <vector>
#include <limits>

#include "from_string.hpp"
#include "Additive.hpp"
#include "VarTuple.hpp"
#include "PrintBlock.hpp"
#include "IntTuple.hpp"
#include "TensorLang.hpp"
#include "interpreter.tab.h"

#define YY_DECL extern "C" int yylex()

%}
%option yylineno
%option nounput
extern int tokval;
%%

[ \t]            ; //ignore whitespace

[\n]            { return EOL; }

#.*  { return COMMENT; }

"p"             { return P; }
"@"             { return ATMARK;  }
"dampening"     { return DAMPENING;  }
"epsilon"       { return EPSILON;  }
"max_iter"      { return MAXITER; }
"engine"        { return ENGINE;  }
"loopy"         { return LOOPY;   }
"brute_force"   { return BRUTEFORCE;  }

-?\.[0-9]+(e-?[0-9]+)?            { yylval.floatPoint = from_string(yytext); return FLOAT; }
(-?[0-9]+)\.([0-9]+)(e-?[0-9]+)?  { yylval.floatPoint = from_string(yytext); return FLOAT; }
-?[0-9]+e-?[0-9]+                 { yylval.floatPoint = from_string(yytext); return FLOAT; }
inf                               { yylval.floatPoint = std::numeric_limits<float>::infinity();  return FLOAT; }
-inf                              { yylval.floatPoint = -std::numeric_limits<float>::infinity();  return FLOAT; }
-[0-9][0-9]*                      { yylval.floatPoint = from_string(yytext); return NEGATIVE_INT; }
[0-9][0-9]*                       { yylval.floatPoint = from_string(yytext); return UNSIGNED_INT; }

\"([^\"])*\"                       { unsigned long len_with_quotes = strlen(yytext); yylval.str = new char[len_with_quotes-2+1]; yylval.str[len_with_quotes-2]=0; memcpy(yylval.str,yytext+1,len_with_quotes-2); return STRING; }

"PMF"           { return PMFTOK;  }
"Pr"            { return PRINT; }
"UNIFORM"       { return UNIFORM; }

"\$" 						{	return ASSIGNTOK; }

"="             { return EQUALS; }
"+"             { return ADD; }
"-"             { return SUB; }
"*"             { return MULT;  }
"/"             { return DIV; }
","             { return COMMA; }
";"             { return SEMICOLON; }

"["       { return LBRACKET; }
"]"       { return RBRACKET; }
"("				{ return LPAREN; }
")"				{ return RPAREN; }

"Gr" { return GRAPH; }

'[^']+'                         { unsigned long len_with_quotes = strlen(yytext); yylval.str = new char[len_with_quotes-2+1]; yylval.str[len_with_quotes-2]=0; memcpy(yylval.str,yytext+1,len_with_quotes-2); return SPECIAL_VARNAME; }

[a-zA-Z_][a-zA-Z0-9_]* { yylval.str = strdup(yytext); return VARNAME; }

.    {std::cerr << "ERROR: syntax error, invalid " << yytext << " on line " << yylineno << std::endl; return INVALID_TOKEN;}
%%
