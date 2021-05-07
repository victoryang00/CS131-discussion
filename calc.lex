/* Mini Calculator */
/* calc.lex */
/* Author: Lan Gao (Modifications by J. Cai, M. Lohstroh, and N. Athreya ) */

%{
#include "heading.h"
#include "calc.tab.h"
void yyerror(const char *s);
static int commentDepth = 0;
%}

digit		[0-9]
int_const	{digit}+

%%

{int_const}	{ if(commentDepth == 0) { yylval.int_val = atoi(yytext); return INTEGER_LITERAL; } }
"+"         { if(commentDepth == 0) {  return PLUS; } }
"-"         { if(commentDepth == 0) {  return MINUS; } }
"*"         { if(commentDepth == 0) { return MULT; } }
"/"         { if(commentDepth == 0) { return DIV; } }
"%"         { if(commentDepth == 0) { return MOD; } }

[ \t]*		{}
[\n]		{ yylineno++;	}

"/\*"       { commentDepth += 1; }
"\*/"       { 
                if (commentDepth > 0) {commentDepth -= 1;} else {  std::cerr << "SCANNER "; yyerror(""); exit(1); }
            }

.           { if (commentDepth == 0) { std::cerr << "SCANNER "; yyerror(""); exit(1); }	}
