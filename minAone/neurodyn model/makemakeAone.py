
import discAone

prob = discAone.Problem
probu = prob.upper()
probl = prob.lower()

nF = discAone.nF

FILE = 'Makefile'

# The name of the IPOPT main file
f = open(FILE,'w')

string = '\\'

f.write('# Makefile for %s IPOPT problem\n\n\
##########################################################################\n\
#    You can modify this example makefile to fit for your own program.   #\n\
#    Usually, you only need to change the five CHANGEME entries below.   #\n\
##########################################################################\n\
\n\
# CHANGEME: This should be the name of your executable\n\
EXE = %s_cpp\n\
\n\
# CHANGEME: Here is the name of all object files corresponding to the source\n\
#           code that you wrote in order to define the problem statement\n\
OBJS = %sminAone_main.o %s\n' % (probl, probl, probl, string))
if nF == 0:
  f.write('       %sminAone_nlp.o\n' % probl)
else:
  f.write('       %sminAone_nlp.o #%s\n' % (probl,string))
  #f.write('       myfunctions.o\n')
f.write('\n\
# CHANGEME: Additional libraries\n\
#ADDLIBS =\n\
\n\
# CHANGEME: Additional flags for compilation (e.g., include flags)\n\
ADDINCFLAGS =\n\
\n\
# CHANGEME: Directory to the sources for the (example) problem definition\n\
# files\n\
#SRCDIR = \n\
#VPATH = \n\
\n\
##########################################################################\n\
#  Usually, you don\'t have to change anything below.  Note that if you   #\n\
#  change certain compiler options, you might have to recompile Ipopt.   #\n\
##########################################################################\n\
\n\
# C++ Compiler command\n\
CXX = clang++\n\
\n\
# C++ Compiler options\n\
CXXFLAGS = -O3 -pipe -DNDEBUG -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long   -DIPOPT_BUILD\n\
\n\
# additional C++ Compiler options for linking\n\
CXXLINKFLAGS = \n\
\n\
# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)\n\
#INCL = `PKG_CONFIG_PATH=/Users/alexanderjulianty/CoinIpopt/build/lib64/pkgconfig:/Users/alexanderjulianty/CoinIpopt/build/lib/pkgconfig:/Users/alexanderjulianty/CoinIpopt/build/share/pkgconfig:  --cflags ipopt` $(ADDINCFLAGS)\n\
INCL = -I`$(CYGPATH_W) /Users/alexanderjulianty/CoinIpopt/build/include/coin`  -I/Users/alexanderjulianty/CoinIpopt/build/include/coin/ThirdParty  -I/Users/alexanderjulianty/CoinIpopt/build/include/coin/ThirdParty   $(ADDINCFLAGS)\n\
\n\
# Directory with header files\n\
IPOPTINCDIR = ${prefix}/include/coin\n\
\n\
# Directory with libipopt.a\n\
IPOPTLIBDIR = ${exec_prefix}/lib\n\
exec_prefix = ${prefix}\n\
prefix = /usr/local/Ipopt\n\
\n\
# Libraries necessary to link with IPOPT\n\
#LIBS = `PKG_CONFIG_PATH=/Users/alexanderjulianty/CoinIpopt/build/lib64/pkgconfig:/Users/alexanderjulianty/CoinIpopt/build/lib/pkgconfig:/Users/alexanderjulianty/CoinIpopt/build/share/pkgconfig:  --libs ipopt`\n\
#LIBS = -link -libpath:`$(CYGPATH_W) /Users/alexanderjulianty/CoinIpopt/build/lib` libipopt.lib  -L/Users/alexanderjulianty/CoinIpopt/build/lib -lcoinhsl -llapack -lblas -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.2.0 -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.2.0/../../.. -lgfortran -lSystem -lquadmath -lm  -L/Users/alexanderjulianty/CoinIpopt/build/lib -lcoinmetis  -framework Accelerate -framework Accelerate -lm  -ldl\n\
LIBS = -L/Users/alexanderjulianty/CoinIpopt/build/lib -lipopt  -L/Users/alexanderjulianty/CoinIpopt/build/lib -lcoinhsl -llapack -lblas -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.2.0 -L/usr/local/Cellar/gcc/7.2.0/lib/gcc/7/gcc/x86_64-apple-darwin16.7.0/7.2.0/../../.. -lgfortran -lSystem -lquadmath -lm  -L/Users/alexanderjulianty/CoinIpopt/build/lib -lcoinmetis  -framework Accelerate -framework Accelerate -lm  -ldl\n\
\n\
# The following is necessary under cygwin, if native compilers are used\n\
CYGPATH_W = echo\n\
\n\
all: $(EXE)\n\
\n\
.SUFFIXES: .cpp .c .o .obj\n\
\n\
$(EXE): $(OBJS)\n\
	bla=;%s\n\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; %s\n\
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(ADDLIBS) $(LIBS)\n\
\n\
clean:\n\
	rm -rf $(EXE) $(OBJS)\n\
\n\
.cpp.o:\n\
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `test -f \'$<\' || echo \'$(SRCDIR)/\'`$<\n\
\n\
\n\
.cpp.obj:\n\
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `if test -f \'$<\'; then $(CYGPATH_W) \'$<\'; else $(CYGPATH_W) \'$(SRCDIR)/$<\'; fi`\n' % (string, string))

f.close()
