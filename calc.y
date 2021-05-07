/* Mini Calculator */
/* calc.y */
/* Author: Lan Gao (Modifications by J. Cai, M. Lohstroh, and N. Athreya ) */

%{
#include "heading.h"
#include <vector>
void yyerror(const char *s);
int yylex(void);
static void print_vec(vector<int>* vec);
%}

%union{
  int		int_val;
  vector<int>* store_vec;
}

%start	input 

%token	<int_val>	INTEGER_LITERAL
%token PLUS MINUS MULT DIV MOD
%type	<int_val>	lit
%type	<int_val>	mult_exp
%type	<int_val>	arith_exp
%type	<int_val>	exp
%type <store_vec> exp_list 

%%

input:	/* empty */
	| exp_list	{ print_vec($1); }
	;

lit: INTEGER_LITERAL { $$ = $1; } ;

mult_exp: lit { $$ = $1; }
        | mult_exp MULT lit { $$ = $1 * $3; }
        | mult_exp DIV lit { $$ = $1 / $3; }
        | mult_exp MOD lit { $$ = $1 % $3; }
        ;

arith_exp: mult_exp { $$ = $1; }
         | arith_exp PLUS mult_exp { $$ = $1 + $3; }
         | arith_exp MINUS mult_exp { $$ = $1 - $3; }
         ;

exp : arith_exp { $$ = $1; } ;

exp_list: exp {$$ = new vector<int>(1,$1);}
	| exp_list exp {$$ = $1; $$->push_back($2);}

%%

void yyerror(const char *s)
{
  extern int yylineno;	// defined and maintained in lex.c
  extern char *yytext;	// defined and maintained in lex.c
  
  cerr << "ERROR: " << string(s) << " at symbol \"" << yytext;
  cerr << "\" on line " << yylineno << endl;
  exit(1);
}

void print_vec(vector<int>* vec) {
  int size = vec->size();
  for (int i = 0; i < size; i++) {
    cout << (*vec)[i] << endl;
  }	
}
