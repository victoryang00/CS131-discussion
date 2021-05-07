/* Mini Calculator */
/* calc.y */
/* Author: Lan Gao (Modifications by J. Cai, M. Lohstroh, and N. Athreya ) */

%{
#include "heading.h"
#include <vector>
void yyerror(const char *s);
int yylex(void);
%}

%union{
  int		int_val;
}

%start	input 
%token	<int_val>	INTEGER_LITERAL
%type	<int_val>	exp
%type	<int_val>	exp_list
%left	PLUS MINUS
%left	MULT DIV MOD

%%

input:	/* empty */
	| exp_list
	;

exp:	INTEGER_LITERAL	{ $$ = $1; }
      | exp PLUS exp	{ $$ = $1 + $3; }
      | exp MULT exp	{ $$ = $1 * $3; }
      | exp MINUS exp	{ $$ = $1 - $3; }
      | exp MOD exp	{ $$ = $1 % $3; }
      | exp DIV exp	{ $$ = $1 / $3; }
	;

exp_list: exp { cout << $1 << endl;}
				| exp_list exp { cout << $2 << endl; }
%%

void yyerror(const char *s)
{
  extern int yylineno;	// defined and maintained in lex.c
  extern char *yytext;	// defined and maintained in lex.c
  
  cerr << "ERROR: " << string(s) << " at symbol \"" << yytext;
  cerr << "\" on line " << yylineno << endl;
  exit(1);
}
