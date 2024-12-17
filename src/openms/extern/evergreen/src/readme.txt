BASIC USE AND EXTENSION OF THE EVERGREEN INFERENCE ENGINE:
==========================

HOW TO #include:
--------------------------
Simply use #include "Evergreen/evergreen.hpp" as done by the demo code. This will include all dependent code. There are one or two optional tools you can #include in addition (as done by the demos).

Demo .cpp files are included in the demos/ directory.

HOW TO COMPILE:
--------------------------
With g++ and clang++: For speed, compile with -O3 -march=native -mtune=native, use minimal error checking #define options (below) and, use fast numeric approximations to pow (below).

The library is header only, and does not need any .lib, .o, or .cpp files to be included.

Do not be surprised if compilation takes some time with aggressive optimizations turned on (e.g., seconds or even a couple of minutes). The TRIOT underpinnings function like a just-in-time (JIT) compiler, which achieves greater speed; however, this requires sometimes more code to be expanded at compile time, which can increasing compilation times.

HOW TO DESIGN YOUR OWN MESSAGE PASSERS:
--------------------------
Simply inherit from either MessagePasser or ContextFreeMessagePasser (see below) and define the necessary abstract functions. You will also need to define an accompanying Dependency type that is associated with your new MessagePasser type. MessagePasser types do the work in the engine, whereas Dependency types provide an interface for creating models (this is not precisely a factory design pattern, but it's similar). 

If your new MessagePasser type has no context (i.e., the edges incident to it do not have distinct meaning from one another), then inherit from ContextFreeMessagePasser. If, on the other hand, the message passer has context (e.g., the pre-multiplication edge and the post-multiplication edge in the ConstantMultiplierMessagePasser), then inherit from MessagePasser.

SCHEDULERS:
==========================
There are two default Scheduler types available: FIFOScheduler and PriorityScheduler (both described below). Both schedulers have parameters defining the dampening; a more aggressive (larger) dampening parameter will force convergence, but may achieve a poorer result. 

FIFOScheduler simply wakes edges when an incident MessagePasser type receives a new message and can pass. Messages are computed lazily (i.e., only once the edge has been selected as the next edge in the ListQueue). Generally, this is lightweight and fast.

PriorityScheduler wakes edges in the same manner, but does not add them to the SetQueue in a lazy manner (i.e., it computes the messages as soon as possible). Since the messages are computed as soon as possible, the PriorityScheduler uses the divergence (via MSE, but this could be generalized) between the old message and the new message, and then in each iteration selects the edge with the highest change (i.e., the least convergent edge). Compared to the FIFOScheduler, this can avoid revisiting edges that are near convergence (but not yet converged) in favor of regions of the graph that are nowhere near convergent. But on the downside, the non-lazy manner of message computation in the PriorityScheduler means that it can be less efficient in the general case. Furthermore, the non-lazy message computation in the PriorityScheduler means that ConvolutionTreeMessagePasser types may not narrow message supports as aggressively in early iterations, because messages are computed as soon as possible, meaning that the message along the final edge incident to a ConvolutionTreeMessagePasser will be computed as soon as messages have been received along the other n-1 edges, and with only n-1 messages received, narrowing cannot be performed. 

INFERENCE GRAPH BUILDERS:
==========================
A particular problem can be composed with several different graphs. One type of graph is the Bethe graph, which creates hyperedges for each variable used and then connects each MessagePasser to the hyperedges corresponding to each variable it uses.

BetheInferenceGraphBuilder is used to construct Bethe graphs, but user-defined inference graph builders can be defined by inheriting from InferenceGraphBuilder and defining the abstract functions. 

ON THE #define OPTIONS:
==========================

These can be declared before the library is #included and will provide error checking at the cost of performance:

ERROR CHECKING: 
--------------------------
#define BOUNDS_CHECK
Verifies every [] operation and lookup (for objects only, it does not check data[index] operations where data is of type T*). It prevents out-of-bounds indexing, but is generally quite expensive.

#define SHAPE_CHECK
Checks shapes once when vectorizing (instead of at every index). Guarantees that all TRIOT expressions will be safe, but much faster than using BOUNDS_CHECK.

#define NUMERIC_CHECK
Checks numeric conditions are met. Fairly inexpensive.

#define ENGINE_CHECK
Checks that parts of the engine itself are working properly. Fairly inexpensive.


SPEEDUPS: 
--------------------------
#define FAST_POW
or
#define FASTER_POW

Replaces some instances of the pow function (e.g., in p-convolution) with a faster numeric approximation. FASTER_POW uses the fastest numeric approximation, but the least accurate. FAST_POW is in between, and using no #define of this category will use the standard pow(...) function. 
