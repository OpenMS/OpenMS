#ifndef _STD_FUNCTIONAL_REPLACEMENT_H
#define _STD_FUNCTIONAL_REPLACEMENT_H
template<class Arg1, class Arg2, class Result>
struct binary_function
{
    using first_argument_type = Arg1;
    using second_argument_type = Arg2;
    using result_type = Result;
};

template <typename ArgumentType, typename ResultType>
struct unary_function
{
    using argument_type = ArgumentType;
    using result_type = ResultType;
};
#endif
