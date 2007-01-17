%{
	// $Id: check_test.l,v 1.2 2006/02/02 18:58:23 marc_sturm Exp $
	// C++ scanner

	#define YY_NO_UNPUT 

	#include <iostream>
	#include <fstream>
	#include <stdio.h>
	#include <stack>
	#include <string>	
	#include <algorithm>
	#include <set>
	#include <list>
	#include <vector>

	using namespace std;
	bool verbose = false;

	string esc(const char* s)
	{
		if (s == 0)
		{
			return "<0x0>";
		}
		string s2(s);
		const char* ptr = s2.c_str();
		string tmp;

		while (*ptr != '\0')
		{
			if (*ptr == '\n')
			{
				tmp += "\\n";
			}
			else if (*ptr == '\t')
			{
				tmp += "\\t";
			}
			else 
			{
				tmp += *ptr;
			}
			++ptr;
		}
		
		return tmp;
	}

	string esc(const std::string& s)
	{
		const char* ptr = s.c_str();
		string tmp;

		while (*ptr != '\0')
		{
			if (*ptr == '\n')
			{
				tmp += "\\n";
			}
			else if (*ptr == '\t')
			{
				tmp += "\\t";
			}
			else 
			{
				tmp += *ptr;
			}
			++ptr;
		}
		
		return tmp;
	}

  string trim(const std::string& s)
  {
  	//search for the begin of truncated string
		std::string::const_iterator begin = s.begin();
		while (begin!=s.end() && (*begin==' ' || *begin=='\t' || *begin=='\n'  || *begin=='\r' ))
		{
			++begin;
		}
		
		//all characters are whitespaces
		if (begin==s.end())
		{
			return "";
		}
		
		//search for the end of truncated string
		std::string::const_iterator end = s.end();
		end--;
		while (end!=begin && (*end==' ' || *end=='\n' || *end=='\t' || *end=='\r' ))
		{
			--end;
		}
		++end;
		
		return string(begin,end);
  }

#define RULE if(debug) cerr << "<" << __LINE__ << "/" << line << "/'" << esc(yytext) << "' -- name= " << name << "  type= " << type << "  args= " << args << "'>" << endl;

#define GOTO(a) \
	if (debug) cerr << "Line " << line << ": <" << __LINE__ << "> -> "#a << "  (name= " << name << "  type= " << type << "  args= " << args << "  yytext= " << yytext << ")" << endl;\
	BEGIN(a)

#define GOSUB(a) \
	if (debug) cerr << "Line " << line << ": gosub -> "#a << endl;\
	return_stack.push(YY_START);\
	BEGIN(a)

#define RETURN \
		if (debug) cerr << "Line " << line <<  ": return" << endl;\
		if (return_stack.empty())\
		{\
			cerr << "Error: return stack empty!" << endl;\
		} else {\
			address=return_stack.top();\
			return_stack.pop();\
		}\
		BEGIN(address)
	// The set continaing the method names found in the class definition.
	std::set<std::string> test_set;

	// The name of the class we are looking for.
	// The first class defined if nothing is specified on the command line.
	std::string test_class = "";

	enum CurrentType
	{
		CLASS,
		NAMESPACE,
		FUNCTION,
		STRUCT,
		ENUM,
		UNION,
		ARRAY,
		TYPEDEF,
		ANY,
		UNKNOWN
	};

	enum Protection
	{
		PUBLIC,
		PRIVATE,
		PROTECTED
	};
	
	static string	name;
	static string	args;
	static string	type;
	static string	program;
	static bool		ignore_namespaces = false;
	static Protection	protection = PUBLIC;
	static int		inner_curly_count = 0;
	static int		skip_curly_start = 0;
	static int		round_bracket_count = 0;
	static int		skip_return   = 0;
	static int		sharp_count   = 0;
	static int		line				 = 1;
	static int		start_line = 0;
	static int		address = 0;
	static bool		debug = false;

	// this stack holds the indentation positions of 
	// opening curly braces
	stack<int>					curly_brackets;
	stack<CurrentType>	current_type;
	stack<string>				current_path;
	stack<int>					return_stack;

	static void print(const char* message, const char* argument = 0)
	{
		if (debug == true)
		{
			cerr << "Line " << line << ": ";
			cerr << message;
			if (argument != 0)
			{
				cerr << " " << argument;
			} 
			cerr << endl;
		}
	}
		
	static void countLine()
	{
		for (const char* c = yytext; *c; ++c)
		{
			line += (*c == '\n');
		}
	}

	static int indentDepth()
	{
		// BAUSTELLE
		return 1;
	}

	static string getPath()
	{
		string path;
		stack<string> tmp = current_path;
		for (; !tmp.empty(); tmp.pop())
		{
			path = tmp.top() + "::" + path;
		}
		
		if (path.size() > 2)
		{
			path.erase(path.size() - 2, 2);
		}
		return path;
	}

	static CurrentType getType()
	{
		if (current_type.empty()) 
			return UNKNOWN;
		else 
			return current_type.top();
	}

	static void pushCurly()
	{
		if (debug) cerr << "Line " << line << ": pushing curly bracket." << endl;
		curly_brackets.push(indentDepth());
	}

	static int popCurly()
	{
		if (debug) cerr << "Line " << line << ": popping curly bracket." << endl;
		int indent = 0;
		if (!curly_brackets.empty())
		{
			indent = curly_brackets.top();
			curly_brackets.pop();
		} 
		else 
		{
			cerr << "Line " << line << ": unable to pop curly_brackets!" << endl;
		}
		return indent;
	}

	static void pushPath(const string& path)
	{
		if (debug == true)	cerr << "Line " << line << ": pushing path: " << path << endl;
		current_path.push(path);
	}

	static void popPath()
	{
		if (!current_path.empty())
		{
			if (debug) cerr << "Line " << line << ": popping path " << current_path.top() << endl;
			current_path.pop();
		} 
		else 
		{
			cerr << "Line " << line << ": unable to pop current_path!" << endl;
		}
	}

	static void pushType(CurrentType type)
	{
		if (debug) cerr << "Line " << line << ": pushing type " << type << endl;
		current_type.push(type);
	}

	static void popType()
	{
		if (!current_type.empty())
		{
			if (debug) cerr << "Line " << line << ": popping type " << current_type.top() << endl;
			current_type.pop();
		} else {
			cerr << "Line " << line << ": unable to pop current_type!" << endl;
		}
	}

	extern "C" int yywrap()
  {
    return 1;
  }

%}

%x	Comment
%x	LineComment
%x	Namespace
%x	UsingNamespace
%x	Define
%x	DefineEnd
%x	DefineEnded
%x	StorageModifier

%x	ClassName
%x  TypeName
%x	Bases

%x	NextSemi
%x	FindMembers
%x  Create
%x	EnumMembers
%x  TypeDef
%x	FindMembersSuffix
%x	GrabSuffixMemo
%x	FindMembersName
%x	Function
%x	Operator
%x	Throws

%x	Array
%x	Round
%x	Curly
%x	SkipCurly
%x	SkipInits
%x	SkipCPP
%x	SkipSemiOrCurly
%x	SkipSemiAndCurly
%x	Sharp
%x	ClassForwardDecl

%%

<*>\x06[^\x06]*\x06	{RULE
	countLine();
	}

"'"[{}]"'"	{RULE
	countLine();
	}

<*>[ \t]*"#i".*	{RULE
		// kill `#ifdef', `#ifndef'
	countLine();
	}

<*>[ \t]*"#e".* {RULE
		// kill `#else', `#endif'
	countLine();
}

<TypeName>[a-zA-Z0-9_::<>&*]+ {RULE
	}

<NextSemi>"{"	{RULE
		// array inialization
		skip_return = NextSemi;
		skip_curly_start = curly_brackets.size();
		GOTO(SkipCurly);
	}

<NextSemi>[;,]	{RULE
		GOTO(FindMembersSuffix);
	}

<EnumMembers>[ \t\n]*[_a-zA-Z][_a-zA-Z0-9]*[ \n\t]*"," {RULE
		if (debug) cerr << "Line " << line << ": definition for enumerator " << yytext << endl;
		countLine();
}

<EnumMembers>[ \t\n]*[_a-zA-Z][_a-zA-Z0-9]*[ \n\t]*=[ \n\t]*[_a-zA-Z][_a-zA-Z0-9:]*[ \n\t]*"," {RULE
		if(debug) cerr << "Line " << line << ": definition for enumerator " << yytext << endl;
		countLine();
}

<EnumMembers>"//".* {RULE
					countLine();
					program += yytext;
					GOSUB(LineComment);
}

<EnumMembers>"[^\n"]*]\" {RULE
					countLine();
					program += yytext;
}

<EnumMembers>"/*"\**[ \t]*			{RULE
					program += yytext;
					GOSUB(Comment);
}

<EnumMembers>"/*"\**[ \t]*\n	{RULE
					program += yytext;
					++line;
					GOSUB(Comment);
}

<EnumMembers>[ \t\n]*"};" {RULE		
		countLine();
		popType();
		popPath();
		name.erase();
		type.erase();
		args.erase();
		GOTO(FindMembers);
}

<EnumMembers>. {RULE
	countLine();
}

<FindMembers>[ \t]*"public"[ \t\n]*":"[ \t\n]* {RULE
		// public section
		protection = PUBLIC;
		type.erase();
		name.erase();
		countLine();
	}

<Create>")"	{RULE
		if ((test_class == "") || (test_class == getPath()))
		{
			// remember the function 
			while ((args.length() > 0)
						 && ((args[args.length() - 1] == ' ') || (args[args.length() -1 ] == '\t')
								 || (args[args.length() - 1] == '\n')))
			{
				args.erase(args.length() - 1, 1);
			}
			test_set.insert(trim(name + args + ")"));
			
		}
		else	
		{
			if (verbose)
			{
				cout << "skipping " << name + args << ") because " << getPath() << " != " << test_class << std::endl;	
			}
		}

		type.erase();
		name.erase();
		args.erase();
		GOTO(FindMembers);
	}

<Create>. {RULE
		args += yytext;
	}
	
<FindMembers>"BALL_CREATE_DEEP"	{RULE
		name = yytext;
		args.erase();
		type.erase();
		GOTO(Create);
	}

<FindMembers>"BALL_CREATE" {RULE
		name = yytext;
		args.erase();
		type.erase();
		GOTO(Create);
	}

<FindMembers>[ \t]*"protected"[ \t\n]*":"[ \t\n]* {RULE
		// protected section
		type.erase();
		name.erase();
		protection = PROTECTED;
		countLine();
	}

<FindMembers>[ \t]*"private"[ \t\n]*":"[ \t\n]*	{RULE
		// private section
		type.erase();
		name.erase();
		protection = PRIVATE;
		countLine();
	}

<FindMembers>[ \t]*"namespace"[ \t\n]+	{RULE
		countLine();
		GOTO(Namespace);
	}

<FindMembers>[ \t]*"using"[ \t\n]+"namespace"[ \t\n]+	{RULE
		// using namespace XXX
		countLine();
		GOTO(UsingNamespace);
	}



<FindMembers>[ \t]*"class"[ \t\n]+[a-z_A-Z0-9]+[ \t]*"<"	{RULE
		// template class forward declaration: ignore
		skip_return = ClassForwardDecl;
		GOTO(Sharp);
	}

<ClassForwardDecl>";" {RULE
		// ignore forward decl
		name.erase();
		type.erase();
		args.erase();
		GOTO(FindMembers);
	}

<FindMembers>[ \t]*"class"[ \t\n]+[a-z_A-Z0-9]+[ \t]*;	{RULE
		// class forward declaration: ignore
		GOTO(FindMembers);
	}

<FindMembers>[ \t]*"class"[ \t\n]+	{RULE
		// class definition
		type.erase();
		name.erase();
		pushType(CLASS);
		countLine();
		GOTO(ClassName);
	}

<FindMembers>typedef[ \n\t]*"struct"[ \t\n]+	{RULE
		// struct definition
		pushType(STRUCT);
		countLine();
		GOTO(ClassName);
	}

<FindMembers>[ \t]*"struct"[ \t\n]+	{RULE
		// struct definition
		pushType(STRUCT);
		countLine();
		GOTO(ClassName);
	}

<FindMembers>[ \t]*"enum"[ \t\n]+"{" {RULE
		// anonymous enum definition
		pushType(ENUM);
		countLine();
		start_line = line;
    pushCurly();
    inner_curly_count = 0;
    pushPath("anonymous enum");
    GOTO(EnumMembers);
	}

<FindMembers>[ \t]*"enum"[ \t\n]+	{RULE
		// enum definition
		pushType(ENUM);
		countLine();
		GOTO(ClassName);
	}

<FindMembers>[ \t]*"typedef"[ \t\n]+[ \t\na-zA-Z0-9<,>:_*&]*";"	{RULE
		// type definition
		if (debug) cerr << "read typedef: " << yytext << endl;
		countLine();
	}

<FindMembers>[ \t]*"union"[ \t\n]+	{RULE
		// union definition
		pushType(UNION);
		countLine();
		GOTO(ClassName);
	}

<FindMembers>[ \t]*"inline"[ \t\n]+	
<FindMembers>[ \t]*"BALL_INLINE"[ \t\n]+	
<FindMembers>[ \t]*"explicit"[ \t\n]+ 
<FindMembers>[ \t]*"auto"[ \t\n]+ 
<FindMembers>[ \t]*"register"[ \t\n]+
<FindMembers>[ \t]*"virtual"[ \t\n]+	
<FindMembers>[ \t]*"extern"[ \t\n]+	
<FindMembers>[ \t]*"mutable"[ \t\n]+	{RULE
		// storage class specifiers
		// and function specifiers
		countLine();
	}

<FindMembers>"operator"/[^a-z_A-Z0-9]	{RULE
		// operator definition
		if (name != "")
		{
			name += " ";
		}
		name += yytext;
		GOTO(Operator);
	}

<Operator>[ \t\n]*"()"?[^(]*		{RULE
		countLine();
		name += yytext;
		GOTO(FindMembersName);
	}

<FindMembers,FindMembersName>[ \t]*":"[ \t]*[0-9]+[ \t]*/";"	{RULE
		// kill obscure bit-width stuff
	}

<FindMembers>::	{RULE
		// append to the current name
		name += yytext;
		GOTO(FindMembersName);
	}

<FindMembers>[a-z_A-Z~.0-9]+	{RULE
		// Normal name
		if (!name.empty())
		{
			if (!type.empty() && ((name[0] != '&') && (name[0] != '*'))) type += " ";
			type += name;
		}
		name = yytext;
	}

<FindMembers>[ \n\t]*"}"[ \n\t]*";" {RULE
		countLine();
		if (getType() != NAMESPACE || !ignore_namespaces)
		{
			popPath();
		}
		popType();
		popCurly();
	}

<FindMembers>"}"	{RULE
		countLine();
		if (getType() == NAMESPACE || !ignore_namespaces)
		{
			popPath();
		}
		popType();
		popCurly();
	}
	

<FindMembersName>"operator"/[^a-z_A-Z0-9] {RULE
		name += yytext;
		GOTO(Operator);
 	}

<FindMembersName>[a-z_A-Z~.0-9]+	{RULE
		name += yytext;
	}

<FindMembersName>.			{RULE
		yyless(0);
		GOTO(FindMembers);
 	}

<FindMembers>^[ \t]*"#"			{RULE
		countLine();
		GOTO(SkipCPP);
	}

<SkipCPP>.
<SkipCPP>"\\"[ \t]*"\n"	{RULE
		countLine();
	}

<SkipCPP>\n				{RULE
		line++;
		GOTO(FindMembers);
	}

<FindMembers>^[ \t]*"#"[ \t]*define[ \t]+ {RULE
		// #define
		GOTO(Define);
	}

<Define>[a-z_A-Z:.0-9]+	{RULE
		name = yytext;
		GOTO(DefineEnd);
	}

<DefineEnd>[ \t]			{RULE
		GOTO(DefineEnd);
	}

<DefineEnd>"("[^)]*")"			{RULE
		countLine();
		args = yytext;
		GOTO(DefineEnded);
	}

<DefineEnd,DefineEnded>"\\"[ \t]*"\n"	{RULE
		line++;
	}

<DefineEnd,DefineEnded>"\n"	{RULE
		line++;
		GOTO(FindMembers);
	}

<FindMembers>[*&]+			{RULE
		name += yytext;
	}

<FindMembers>[;=,]			{RULE
		GOTO(FindMembersSuffix);
		if ((name != "") && (type != ""))
		{
			if (debug) cerr << "Line "<< line << ": member " << type << " " << getPath() << name << endl;
		} else {
			popPath();
			popType();
		}

		name.erase();
		args.erase();
		type.erase();

		if(*yytext == '=') 
		{
			GOTO(NextSemi);
		}
	}

<FindMembers>"["			{RULE
		args += yytext;
		sharp_count = 1;
		GOTO(Array);
	}

<Array>"]"	{RULE
		args += *yytext;
		if(--sharp_count <= 0)	GOTO(FindMembers);
	}

<Array>"["	{RULE
		args += *yytext;
		sharp_count++;
	}

<Array>.				{RULE
		countLine();
		args += *yytext;
	}

<FindMembers>"<"	{RULE
		if (!type.empty() && ((name[0] != '&') && (name[0] != '*'))) type += " ";
		type += name;
		name.erase();
		sharp_count++;
		skip_return = FindMembers;
		type += "<";
		GOTO(Sharp);
	}

<Sharp>">"	{RULE
		type += yytext;
		if(--sharp_count <= 0) GOTO(skip_return);
	}

<Sharp>"<"	{RULE
					type += yytext;
					sharp_count++;
					}

<Sharp>.	{RULE
					countLine();
					type += *yytext;
					}

<Curly>[^\n{}"/]*			{RULE
					countLine();
					program += yytext;
					}

<Curly>"//".*				{RULE
		countLine();
					program += yytext;
					GOSUB(LineComment);
					}

<Curly>"[^\n"]*]\"			{RULE
					countLine();
					program += yytext;
					}

<Curly>"/*"\**[ \t]*			{RULE
					program += yytext;
					GOSUB(Comment);
					}

<Curly>"/*"\**[ \t]*\n	{RULE
					program += yytext;
					++line;
					GOSUB(Comment);
					}

<Curly>"{"	{RULE
		program += yytext;
		pushCurly();
		++inner_curly_count;
		pushType(ANY);
		pushPath("ANY");
	}

<Curly>"}"	{RULE
		popCurly();
		popPath();
		popType();
		countLine();
		if (inner_curly_count == 0)	
		{
			type.erase();
			args.erase();
			GOTO(NextSemi);
		}
	}		
				

<Curly>\n	{RULE
		program += yytext;
		line++;
	}

<Curly>.	{RULE
		countLine();
		program += yytext;
	}

<FindMembers>"("			{RULE
		args = yytext;
		GOSUB(Round);
	}

<FindMembers>")"	{RULE
		// back from Round
		GOTO(Function);
	}

<Round>"("	{RULE
		args += *yytext;
		++round_bracket_count;
	}

<Round>")" {RULE
		args += *yytext;
		if(round_bracket_count)
		{
			--round_bracket_count;
		}	
		else
		{
			yyless(0);
			RETURN;
		}
	}


<Round>[ \t\n]*","[ \t\n]*	{RULE
		countLine();
		args += ", ";
	}

<Round>[ \t\n]+	{RULE
		countLine();
		args += ' ';
	}

<Round>.				{RULE
		countLine();
		args += *yytext;
	}

<Function>[ \n\t]*"const"[ \t\n]*	{RULE
					args += " const ";
					countLine();
					}

<Function>[ \n\t]*"throw"[ \t\n]*"("	{RULE
		args += " throw(";
		countLine();
		GOSUB(Round);
	}

<Function>[ \n\t]*"const"[ \t\n]*"throw"[ \t\n]*"("	{RULE
		args += " const throw(";
		countLine();
		GOSUB(Round);
	}

<Function>"("	{RULE
		type += name;
		name = args;
		args = yytext;
		GOSUB(Round);
	}

<Function>[":"";""{"]			{RULE
		// document public methods only
		if ((protection == PUBLIC) && (getType() == CLASS))
		{
			if ((test_class == "") || (test_class == getPath()))
			{
				// remember the function 
				while ((args.length() > 0)
							 && ((args[args.length() - 1] == ' ') || (args[args.length() -1 ] == '\t')
									 || (args[args.length() - 1] == '\n')))
				{
					args.erase(args.length() - 1, 1);
				}
									 
				if (!type.empty()) type += " ";
				test_set.insert(trim(type + name + args));
			}
			else	
			{
				if (verbose)
				{
					cout << "skipping " << name + args << " because " << getPath() << " != " << test_class << std::endl;
				}
			}
		}

		name.erase();
		type.erase();
		args.erase();
		if(*yytext == '{')
		{
			pushCurly();
			skip_return = FindMembers;
			skip_curly_start = curly_brackets.size();
			GOTO(SkipCurly);
		}
		else
		{
			if(*yytext == ':')
			{
				GOTO(SkipInits);
			}
			else
			{
				GOTO(FindMembersSuffix);
			}
		}
	}

<Function>.				{RULE
	countLine();
					}

<FindMembersSuffix>\n			{RULE
		yyless(0);
		GOTO(FindMembers);
	}

<FindMembersSuffix>.*"//"	{RULE
		countLine();
		yyless(0);
		GOTO(FindMembers);
	}

<FindMembersSuffix>.*"/\*"	{RULE
		countLine();
		yyless(0);
		GOTO(FindMembers);
	}

<FindMembersSuffix>.	{RULE
		yyless(0);
		GOTO(FindMembers);
	}
				     
<SkipSemiOrCurly>"{"	{RULE
		skip_return = FindMembers;
		skip_curly_start = curly_brackets.size();
		GOTO(SkipCurly);
	}

<SkipSemiOrCurly>;	{RULE
		yyless(0);
		GOTO(FindMembers);
	}

<SkipCurly>"{"				{RULE
		// skip this one...
		pushCurly();
	}

<SkipCurly>"}"				{RULE
		if ((int)curly_brackets.size() > skip_curly_start)
		{
			popCurly();
		}
		else
		{
			GOTO(skip_return);
		}
	}

<SkipInits>"{"	{RULE
		skip_return = SkipSemiAndCurly;
		skip_curly_start = curly_brackets.size();
		GOTO(SkipCurly);
	}

<SkipInits>";"	{RULE
		GOTO(FindMembers);
	}

<SkipSemiAndCurly>[ \t\n]	{RULE
					yyless(0);
					GOTO(FindMembers);
	}

<SkipSemiAndCurly>.			{RULE
		// i.e. a semicolon (or anything else)
					if(*yytext != ';')
					    yyless(0);
					countLine();
					GOTO(FindMembers);
					}

<Bases>";"				
	{
		GOTO(FindMembersSuffix);
	}

<ClassName>[a-z_A-Z0-9\[\]*&]+		{RULE
		if ((name.length() > 0) && !current_type.empty() && (current_type.top() == UNION))
		{
			print("found union" , yytext);
		}
		else if ((name.length() > 0) && !current_type.empty() && (current_type.top() == ENUM))
		{
			print("found enum" , yytext);
		}
		else if ((name.length() > 0) && !current_type.empty() && (current_type.top() == STRUCT))
		{
			print("found struct" , yytext);
		}
		if (!type.empty())
		{
			type += " ";
		}
		type += name;
		name = yytext;
		pushPath(name);
		if (verbose)
		{
			cout << "reading entries for class " << getPath() << endl;
		}
	}

<ClassName>"<" {RULE
		sharp_count++;
		skip_return = ClassName;
		GOTO(Sharp);
	}

<ClassName>[ \t]*":"[ \t]*	{RULE
		args = ":";
		GOTO(Bases);
	}

<Bases,ClassName>[ \t]*"{"[ \t]*	{RULE
		start_line = line;
		pushCurly();
		inner_curly_count = 0;
		if (getType() == ENUM)
		{
			GOTO(EnumMembers);
		} 
		else 
		{
			GOTO(FindMembers);
		}
	}

<Bases,ClassName>[ \t]*";"	{RULE
		start_line = line;
		if (debug) cerr << "Line " << line << ": found forward class declaration of " << getPath() << yytext << endl;
		inner_curly_count = 0;
		GOTO(FindMembers);
	}

<Bases>[a-z_A-Z*.<>0-9:]+		{RULE
					countLine();
					args += ' ';
					args += yytext;
					}

<Bases>","				{RULE
					args += ',';
					}

<LineComment>[^\n]*\n {RULE
		countLine();
		RETURN;
	}

<Comment>\n				{RULE
					program += yytext;
					line++;
					}

<Comment>.				{RULE
					countLine();
					program += yytext;
					}

<Comment>.*"*/"				{RULE
					countLine();
					program += yytext;
					RETURN;
					}

<FindMembers>[ \t\n]*"/\*" {RULE
					countLine();
					GOSUB(Comment);
					}

<FindMembers>[ \t\n]*"//" {RULE
					countLine();
					GOSUB(LineComment);
					}

<Namespace>[a-z_A-Z0-9]*		{RULE
		name = yytext;
		print("found namespace ", yytext);
		pushType(NAMESPACE);
		if (!ignore_namespaces)
		{
			pushPath(name);
		}
		type.erase();
		name.erase();
		args.erase();
	}

<Namespace>"{"				{RULE
		pushCurly();
		GOTO(FindMembers);
	}

<UsingNamespace>[a-z_A-Z0-9"::"]*	{RULE
		// using namespace : ignored
	}

<UsingNamespace>;	{RULE
		GOTO(FindMembers);
	}

<*>"\n"	{RULE
		countLine();
	}

<*>"}" {RULE
		countLine();
	}

<*>. {RULE
		countLine();
	}


%%

set<string> tests_found;
void parseTestFile(const char* filename)
{
	ifstream is(filename);
	string line;
	while (is)
	{
		getline(is, line);
		string::size_type check_pos = line.find("CHECK(");
		if (check_pos != string::npos)
		{
			line.erase(0,check_pos+6);
			line.erase(line.length() - 1, 1);

			// Remove round brackets
			while ((line.size() > 0) && (line[0] == '(') && (line[line.size() - 1] == ')'))
			{
				line.erase(0, 1);
				line.erase(line.length() - 1, 1);
			}

			if (tests_found.count(line) > 0)
			{
				cout << "WARNING: duplicate test for " << line << " of " << filename << endl;
			}
			tests_found.insert(trim(line));
		}
	}
	is.close();
}

int main(int argc, char** argv)
{
	FILE* infile;
	yyin = stdin;
	yyout = stdout;
	
	if (argc < 2)
	{	
		cerr << "check_test <header> [<test_file> [<class>]] [-a] [-c] [-d]\n\n"
				 << "  Parse the header file (and optionally a test file) and\n"
				 << "  determine which tests are required for a specific class.\n\n"
				 << "    -a : print all methods defined in the class.\n"
         << "         The default is to print the missing ones only if <test_file>\n"
				 << "         is given.\n\n"
				 << "    -c : create output for the test file directly.\n"
         << "         This option causes the list of missing methods to be printed\n"
         << "         in the correct format (CHECK()/RESULT) for the test file.\n\n"
         << "    -v : verbose mode -- print some information on what's going on.\n\n"
         << "    -d : debug mode -- if you want to debug the parser, be my guest!\n\n"
         << "  <test_file> is a BALL/BioMiner/OpenMS-compatible test file.\n"
         << "         check_test searches for lines beginning with 'CHECK(' only.\n"
         << "         The method names have to match *exactly* (i.e. char by char)\n"
				 << "         in order to be recognized. Additional tests are ignored.\n"
         << "         The order of the tests is irrelevant.\n\n"
         << "  <class> is the name of a specific class. The methods of all other classes\n"
				 << "         declared in the header file as well are ignored. The class name has\n"
         << "         has to be fully qualified if it is inside a namespace or if it is a\n"
         << "         nested class.\n" 
				 << endl;

    return 3; 																																					
	}

	infile = fopen(argv[1], "r");
	if (infile)
	{
		yyin = infile;
	}
	else
	{
		cerr << "Cannot open file " << argv[1] << "!" << endl;
		return 4;
	}

	bool create_test_bodies = false;
	bool print_all_methods = false;

	// We have a class specified as well.
	for (int i = 2; i < argc; ++i)
	{
		if (!strcmp(argv[i], "-c"))
		{
			create_test_bodies = true;
		}
		else if (!strcmp(argv[i], "-d"))
		{
			debug = true;
		}
		else if (!strcmp(argv[i], "-v"))
		{
			verbose = true;
		}
		else if (!strcmp(argv[i], "-a"))
		{
			print_all_methods = true;
		}
		else if ((argv[i][0] != '-') && (i > 2))
		{
			test_class = argv[i];
		}
		else if ((i > 2) && (argv[i][0] == '-'))
		{
			cerr << "unknown command line option: " << argv[i] << endl;
			return 1;
		}
	}
		
	// Parse the header file
	BEGIN(FindMembers);
	yylex();

	// Parse the test file.
	if ((argc > 2) && (argv[2] != "-c"))	
	{
		parseTestFile(argv[2]);
	}
	
	set<string> missing;
	if (print_all_methods)
	{
		cout << "Methods required for test of " << argv[1] << endl;
		copy(test_set.begin(), test_set.end(), inserter(missing, missing.begin()));
	}
	else
	{
		set_difference(test_set.begin(), test_set.end(),
									 tests_found.begin(), tests_found.end(),
									 inserter(missing, missing.begin()));
	}

	set<string>::iterator it(missing.begin());
	for (; it != missing.end(); ++it)
	{
		if (create_test_bodies)
		{
			if (it->find(',') == string::npos)
			{		
				cout << "CHECK(" << *it << ")\n"
						 << "  // ???\n"
						 << "RESULT\n"
						 << endl;
			}
			else	
			{
				cout << "CHECK((" << *it << "))\n"
						 << "  // ???\n"
						 << "RESULT\n"
						 << endl;
			}
		}
		else	
		{
			cout << "  '" << esc(*it) << "'" << endl;
		}
	}

  if (!print_all_methods)
  {
    missing.clear();
    set_difference(tests_found.begin(), tests_found.end(),
                   test_set.begin(), test_set.end(),
                   inserter(missing, missing.begin()));

    if (!missing.empty() && (argc > 2))
    {
			//check if all missing methods are [EXTRA]
			vector<string> no_extra;
			for (set<string>::iterator it = missing.begin(); it != missing.end(); ++it)
      {
				if (it->find("[EXTRA]") == string::npos)
				{
					no_extra.push_back(*it);
				}
			}
			if (no_extra.size()>0)
			{
	      cout << "Unknown methods in " << argv[2] << endl;
	      for (vector<string>::iterator it = no_extra.begin(); it != no_extra.end(); ++it)
	      {
					cout << "- '" << esc(*it) << "'" << endl;
				}
			}
		}
	}

	return 0;
}
