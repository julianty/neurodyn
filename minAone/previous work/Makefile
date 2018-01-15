# Makefile for simple_nakl IPOPT problem

##########################################################################
#    You can modify this example makefile to fit for your own program.   #
#    Usually, you only need to change the five CHANGEME entries below.   #
##########################################################################

# CHANGEME: This should be the name of your executable
EXE = simple_nakl_cpp

# CHANGEME: Here is the name of all object files corresponding to the source
#           code that you wrote in order to define the problem statement
OBJS = simple_naklminAone_main.o \
       simple_naklminAone_nlp.o

# CHANGEME: Additional libraries
#ADDLIBS =

# CHANGEME: Additional flags for compilation (e.g., include flags)
ADDINCFLAGS =

# CHANGEME: Directory to the sources for the (example) problem definition
# files
#SRCDIR = 
#VPATH = 

##########################################################################
#  Usually, you don't have to change anything below.  Note that if you   #
#  change certain compiler options, you might have to recompile Ipopt.   #
##########################################################################

# C++ Compiler command
CXX = clang++

# C++ Compiler options
CXXFLAGS = -O3 -pipe -DNDEBUG -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long   -DIPOPT_BUILD

# additional C++ Compiler options for linking
CXXLINKFLAGS = 

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
#INCL = `PKG_CONFIG_PATH=/Users/alexanderjulianty/CoinIpopt/build/lib64/pkgconfig:/Users/alexanderjulianty/CoinIpopt/build/lib/pkgconfig:/Users/alexanderjulianty/CoinIpopt/build/share/pkgconfig:  --cflags ipopt` $(ADDINCFLAGS)
INCL = -I`$(CYGPATH_W) /Users/alexanderjulianty/CoinIpopt/build/include/coin`  -I/Users/alexanderjulianty/CoinIpopt/build/include/coin/ThirdParty  -I/Users/alexanderjulianty/CoinIpopt/build/include/coin/ThirdParty   $(ADDINCFLAGS)

# Directory with header files
IPOPTINCDIR = ${prefix}/include/coin

# Directory with libipopt.a
IPOPTLIBDIR = ${exec_prefix}/lib
exec_prefix = ${prefix}
prefix = /usr/local/Ipopt

# Libraries necessary to link with IPOPT
#LIBS = `PKG_CONFIG_PATH=/Users/alexanderjulianty/CoinIpopt/build/lib64/pkgconfig:/Users/alexanderjulianty/CoinIpopt/build/lib/pkgconfig:/Users/alexanderjulianty/CoinIpopt/build/share/pkgconfig:  --libs ipopt`
#LIBS = -link -libpath:`$(CYGPATH_W) /Users/alexanderjulianty/CoinIpopt/build/lib` libipopt.lib  -L/Users/alexanderjulianty/CoinIpopt/build/lib -lcoinhsl -llapack -lblas -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.2.0 -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.2.0/../../.. -lgfortran -lSystem -lquadmath -lm  -L/Users/alexanderjulianty/CoinIpopt/build/lib -lcoinmetis  -framework Accelerate -framework Accelerate -lm  -ldl
LIBS = -L/Users/alexanderjulianty/CoinIpopt/build/lib -lipopt  -L/Users/alexanderjulianty/CoinIpopt/build/lib -lcoinhsl -llapack -lblas -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.2.0 -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.2.0/../../.. -lgfortran -lSystem -lquadmath -lm  -L/Users/alexanderjulianty/CoinIpopt/build/lib -lcoinmetis  -framework Accelerate -framework Accelerate -lm  -ldl

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(ADDLIBS) $(LIBS)

clean:
	rm -rf $(EXE) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `test -f '$<' || echo '$(SRCDIR)/'`$<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `if test -f '$<'; then $(CYGPATH_W) '$<'; else $(CYGPATH_W) '$(SRCDIR)/$<'; fi`
