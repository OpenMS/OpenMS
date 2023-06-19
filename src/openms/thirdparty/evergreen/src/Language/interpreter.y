%{
#include <iostream>
#include <unordered_map>
#include <string>
#include <cstring>
#include <algorithm>
#include <vector>
#include <limits>

double p = 1;

#include "Additive.hpp"
#include "VarTuple.hpp"
#include "PrintBlock.hpp"
#include "IntTuple.hpp"
#include "TensorLang.hpp"
#include "../Evergreen/evergreen.hpp"
#include "LangEngine.hpp" // fixme: rename

#define YYDEBUG 1
#define YYMAXDEPTH  100000
#define YYINITDEPTH 100000

  const double default_dampening = 0;
  const double default_epsilon = 1e-15;
  const long default_max_iter = 1L<<14;
  LangEngine engine(default_dampening, default_epsilon, default_max_iter);

  extern "C" FILE *yyin;
  extern "C" int yylex();
  extern "C" int yyparse();
  extern int yylineno;
  int N_index = 0;
  void yyerror(const char *s);
%}

%union {
  double floatPoint;
  int integer;
  char*str;
  Additive*additive;
  VarTuple*varTuple;
  IntTuple*intTuple;
  TensorLang*tensor;
  PrintBlock*printBlock;
}

%error-verbose

%token EOL
%token COMMENT
%token PMFTOK
%token UNIFORM
%token PRINT
%token GRAPH
%token SEMICOLON
%token LPAREN RPAREN
%token LBRACKET RBRACKET

%token VARNAME
%token SPECIAL_VARNAME

%token NEGATIVE_INT
%token UNSIGNED_INT
%token FLOAT
%token ASSIGNTOK
%token ADD SUB MULT DIV EQUALS COMMA

%token STRING

%token ENGINE
%token DAMPENING
%token EPSILON
%token MAXITER
%token LOOPY
%token BRUTEFORCE

%token P
%token ATMARK

%token INVALID_TOKEN

%type<floatPoint> NEGATIVE_INT
%type<floatPoint> UNSIGNED_INT
%type<floatPoint> FLOAT
%type<floatPoint> floatInput
%type<intTuple> intTupleExpr

%type<integer> intExpr
%type<floatPoint> floatExpr
%type<floatPoint> floatAddSubExpr
%type<floatPoint> floatMultDivExpr
%type<floatPoint> atomicExpr

%type<str> STRING

%type<str> VAR
%type<str> VARNAME
%type<str> SPECIAL_VARNAME

%type<varTuple> varTuple

%type<str> pmfLine
%type<str> uniformPMF
%type<str> categoricalPMF

%type<additive> addSubExpr
%type<str> additiveLine

%type<str> printLine
%type<printBlock> printBlock

%type<tensor> floatTensor
%type<tensor> floatTensorTuple
%type<tensor> floatTupleExpr
%%

prog: multiLine
;

multiLine: multiLine EOL line
| line
;

line: additiveLine
| pmfLine 
| printLine
| saveGraphLine
| commentLine
| pLine
| engineLine
| error
| INVALID_TOKEN
| /* NULL */
;

engineLine: ATMARK ENGINE EQUALS engineType
;

engineType: BRUTEFORCE LPAREN RPAREN {
  engine.set_engine(new BruteForceInferenceEnginesBuilder);
}

| LOOPY LPAREN RPAREN {
  engine.set_engine(new BeliefPropagationInferenceEnginesBuilder(default_dampening, default_epsilon, default_max_iter));
}
| LOOPY LPAREN ATMARK DAMPENING EQUALS floatExpr COMMA ATMARK EPSILON EQUALS floatExpr COMMA ATMARK MAXITER EQUALS intExpr RPAREN {
  double dampening = $6;
  assert(dampening >= 0.0 && dampening <= 1.0);
  double epsilon = $11;
  assert(epsilon >= 0.0 && epsilon <= 1.0);
  long max_iter = $16;
  assert(max_iter >= 0.0);

  engine.set_engine(new BeliefPropagationInferenceEnginesBuilder(dampening, epsilon, max_iter));
}
;
 
pLine: ATMARK P EQUALS floatInput {
  double p_value = $4;
  assert(p_value>0);
  p = $4;
}
;

commentLine: line COMMENT
;

additiveLine: addSubExpr EQUALS addSubExpr {
  std::vector<std::vector<std::string> > left_vars;
  std::vector<std::vector<std::string> > right_vars;
  for (std::vector<std::string> plus_var : $1->plus_vars)
    left_vars.push_back(plus_var);
  for (std::vector<std::string> minus_var : $1->minus_vars)
    right_vars.push_back(minus_var);
  for (std::vector<std::string> plus_var : $3->plus_vars)
    right_vars.push_back(plus_var);
  for (std::vector<std::string> minus_var : $3->minus_vars)
    left_vars.push_back(minus_var);

  unsigned long additive_length = left_vars[0].size();
  for (const std::vector<std::string> & tup : left_vars)
    if (tup.size() != additive_length)
      std::cerr << "ERROR: additive dependency error, Not all tuples of the additive share the same length" << std::endl;
  for (const std::vector<std::string> & tup : right_vars)
    if (tup.size() != additive_length)
      std::cerr << "ERROR: additive dependency error, Not all tuples of the additive share the same length" << std::endl;

  if(left_vars.size() == 1)
    engine.insert_dependency(new AdditiveDependency<std::string>(right_vars,left_vars[0],p));
  else if(right_vars.size()==1)
    engine.insert_dependency(new AdditiveDependency<std::string>(left_vars,right_vars[0],p));
  else {
    engine.insert_dependency(new AdditiveDependency<std::string>(left_vars,{"@N" + to_string(N_index)},p));
    engine.insert_dependency(new AdditiveDependency<std::string>(right_vars,{"@N" + to_string(N_index)},p));
    ++N_index;
  }
  // we're finished with both addSubExpr types:
  delete $1;
  delete $3;
}
| addSubExpr EQUALS addSubExpr ATMARK P EQUALS floatExpr {
  double p_value = $7;
  assert(p_value>0);
  std::vector<std::vector<std::string> > left_vars;
  std::vector<std::vector<std::string> > right_vars;
  for (std::vector<std::string> plus_var : $1->plus_vars)
    left_vars.push_back(plus_var);
  for (std::vector<std::string> minus_var : $1->minus_vars)
    right_vars.push_back(minus_var);
  for (std::vector<std::string> plus_var : $3->plus_vars)
    right_vars.push_back(plus_var);
  for (std::vector<std::string> minus_var : $3->minus_vars)
    left_vars.push_back(minus_var);

  unsigned long additive_length = left_vars[0].size();
  for (const std::vector<std::string> & tup : left_vars)
    if (tup.size() != additive_length)
      std::cerr << "ERROR: additive dependency error, Not all tuples of the additive share the same length" << std::endl;
  for (const std::vector<std::string> & tup : right_vars)
    if (tup.size() != additive_length)
      std::cerr << "ERROR: additive dependency error, Not all tuples of the additive share the same length" << std::endl;

  if(left_vars.size() == 1)
    engine.insert_dependency(new AdditiveDependency<std::string>(right_vars,left_vars[0],$7));
  else if(right_vars.size()==1)
    engine.insert_dependency(new AdditiveDependency<std::string>(left_vars,right_vars[0],$7));
  else {
    engine.insert_dependency(new AdditiveDependency<std::string>(left_vars,{"@N1"},$7));
    engine.insert_dependency(new AdditiveDependency<std::string>(right_vars,{"@N2"},$7));
  }
  // we're finished with both addSubExpr types:
  delete $1;
  delete $3;
}
;

pmfLine: uniformPMF
| categoricalPMF
;

uniformPMF: PMFTOK LPAREN varTuple RPAREN LPAREN intTupleExpr RPAREN LPAREN intTupleExpr RPAREN UNIFORM {
  std::vector<std::string> labels($3->vars.begin(), $3->vars.end());
  if ($6->ints.size() != $9->ints.size() || $6->ints.size() != labels.size() || labels.size() != $9->ints.size())
    std::cerr << "ERROR: PMF error, number of variables, low bounds and high bounds do not match on line " << yylineno << std::endl;
  std::vector<long> _shape;
  long flat_length = 1;
  for (unsigned long i=0; i<$6->ints.size(); ++i)
    _shape.push_back($9->ints[i] - $6->ints[i] + 1);
  for (long dim : _shape)
    flat_length *= dim; 
  std::vector<double> _flat_vect(flat_length, 1.0/double(flat_length));
  Vector<double> flat_vect(_flat_vect);
  Vector<long> shape(_shape);
  Tensor<double> tensor(shape,flat_vect);

  PMF pmf($6->ints,tensor);
  LabeledPMF<std::string> lpmf(labels,pmf);
  TableDependency<std::string> table(lpmf, p);
  engine.insert_dependency(new TableDependency<std::string>(table));

  delete $6;
  delete $9;
  delete $3;
}
;

categoricalPMF: PMFTOK LPAREN varTuple RPAREN LPAREN intTupleExpr RPAREN floatTensor {
  std::vector<std::string> labels($3->vars.begin(), $3->vars.end());  
  Vector<double> flat_vect($8->flat_vector);
  std::reverse($8->shape.begin(),$8->shape.end());
  Vector<long> shape($8->shape);

  Tensor<double> tensor(shape,flat_vect);

  if (labels.size() != tensor.dimension() || tensor.dimension() != $6->ints.size() || labels.size() != $6->ints.size())
    std::cerr << "ERROR: PMF error, number of variables, offsets, and dimension of tensor do not match on line " << yylineno << std::endl;

  PMF pmf($6->ints,tensor);  
  LabeledPMF<std::string> lpmf(labels,pmf);
  TableDependency<std::string> table(lpmf, p);
  engine.insert_dependency(new TableDependency<std::string>(table));

  delete $8;
  delete $6;
  delete $3;
}
| PMFTOK LPAREN varTuple RPAREN LPAREN intTupleExpr RPAREN floatTensor ATMARK P EQUALS floatExpr {
  double p_value = $12;
  assert(p_value>0);
  Vector<double> flat_vect($8->flat_vector);
  std::reverse($8->shape.begin(),$8->shape.end());
  Vector<long> shape($8->shape);
  Tensor<double> tensor(shape,flat_vect);
  PMF pmf($6->ints,tensor);
  std::vector<std::string> labels($3->vars.begin(), $3->vars.end());
  LabeledPMF<std::string> lpmf(labels,pmf);
  TableDependency<std::string> table(lpmf,$12);
  engine.insert_dependency(new TableDependency<std::string>(table));

  delete $6;
  delete $3;
}
;

printLine: PRINT LPAREN printBlock RPAREN {
  std::vector<std::vector<std::string> > result_vars = $3->vars;
  engine.print($3->vars);
  // we are finished with the print block now:
  delete $3;
}
| PRINT LPAREN STRING RPAREN {
  printf("%s\n", $3);
  delete $3;
}
| PRINT LPAREN RPAREN {
  engine.print_normalization_constant();
}
| PRINT LPAREN MULT RPAREN {
  engine.recompute_and_print_normalization_constant();
}
;

printBlock: printBlock SEMICOLON varTuple {
  $$->vars.push_back($3->vars);
  
  // we're now finished with the varTuple:
  delete $3;
}
| varTuple {
  $$ = new PrintBlock;
  $$->vars.push_back($1->vars);
  delete $1;
} 
;

saveGraphLine: GRAPH LPAREN STRING RPAREN {
  engine.save_graph($3);
  // todo: is there a better way to do this without an ugly system call?
  std::string instruction_a = std::string("python $EVERGREEN_SRC_PATH/Utility/draw_dot.py ") + $3 + std::string(" ") + $3 + std::string(".png");
  system(instruction_a.c_str());
  // todo: hard-coded for linux... need python executable that draws
  // the graph in interactive mode
  std::string instruction_b = std::string("xdg-open ") + $3 + std::string(".png");
  system(instruction_b.c_str());
}

floatTensor: LBRACKET floatTensorTuple RBRACKET {
  $$ = $2;
}
| LBRACKET floatTupleExpr RBRACKET {
  $$ = $2;
}
;

floatTensorTuple: floatTensorTuple COMMA floatTensor {
  long array_in_tensor_length = $1->shape[$$->shape.size()-2];
  long array_2_in_tensor_length = $3->shape[$$->shape.size()-2];
  assert(array_in_tensor_length == array_2_in_tensor_length);
  $$ = $1;  
  for (double flt: $3->flat_vector)
    $$->flat_vector.push_back(flt);
  $$->shape[$$->shape.size()-1] = $$->shape[$$->shape.size()-1]+1;
  delete $3;
}
| floatTensor {
  $$->shape.push_back(1);
}
;

floatTupleExpr: floatTupleExpr COMMA floatExpr {
  $$->flat_vector.push_back($3);
  $$->shape[0] = $$->shape[0]+1;
}
| floatExpr {
  $$ = new TensorLang;
  $$->shape = {1};
  $$->flat_vector.push_back($1);
}
;

addSubExpr: addSubExpr ADD varTuple {
  $$->plus_vars.push_back($3->vars);
  delete $3;
}
| addSubExpr SUB varTuple {
  $$->minus_vars.push_back($3->vars);
  delete $3;
}
| varTuple {
  $$ = new Additive;
  $$->plus_vars.push_back($1->vars);
  delete $1;
}
;

intTupleExpr: intTupleExpr COMMA intExpr {
  $$ = $1;
  $$->ints.push_back($3);
}
| intExpr {
  $$ = new IntTuple;
  $$->ints.push_back($1);
}
;

varTuple: varTuple COMMA VAR {
  $$ = $1;
  $$->vars.push_back($3);
  free($3);
}
| VAR {
  $$ = new VarTuple;
  $$->vars.push_back($1);
  free($1);
}
;

VAR: VARNAME
| SPECIAL_VARNAME

/********************************
Expressions and Numbers:
********************************/

intExpr: floatExpr{
  $$ = (long)$1;
}
;

floatExpr: floatAddSubExpr
;

floatAddSubExpr: floatAddSubExpr ADD floatMultDivExpr {
  $$ = $1 + $3;
}
| floatAddSubExpr SUB floatMultDivExpr {
  $$ = $1 - $3;
}
| floatMultDivExpr
;

floatMultDivExpr: floatMultDivExpr MULT atomicExpr {
	$$ = $1 * $3;
}
| floatMultDivExpr DIV atomicExpr {
	if ($3 == 0) {
	  std::cerr << "ERROR: Cannot divide by 0!" << std::endl;
	  exit(1);
	}
	else {
  	  $$ = $1 / $3;
	}
}
| atomicExpr
;

atomicExpr: floatInput
| LPAREN floatExpr RPAREN {
	$$ = $2;
}
;

floatInput: FLOAT
| NEGATIVE_INT
| UNSIGNED_INT
;

%%

int main(int argc, char **argv) {
  FILE *fh;
  if (argc == 2 && (fh = fopen(argv[1], "r")))
    yyin = fh;
  //yydebug = 1; (uncomment to enable debugging)
  yyparse();

  //  fclose(fh);
  return 0;
}

/* Display error messages */
void yyerror(const char *s) {
  
  fprintf(stderr, "ERROR: %s on line %d\n", s, yylineno);
}
