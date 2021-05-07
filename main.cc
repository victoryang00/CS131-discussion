/* main.cc */

#include "heading.h"
#include <cstdlib>
#include <vector>

// prototype of bison-generated parser function
int yyparse();
std::vector <int> output_list;

int main(int argc, char **argv)
{
  if ((argc > 1) && (freopen(argv[1], "r", stdin) == NULL))
  {
    cerr << argv[0] << ": File " << argv[1] << " cannot be opened.\n";
    exit( 1 );
  }
  
  yyparse();

  return 0;
}


