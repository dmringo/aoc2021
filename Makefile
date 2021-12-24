 
ifdef CLANG_BUILD
CXX = clang++
CXX_FLAGS = -Wall -std=c++17 -o $@.exe
ifdef DEBUG 
CXX_FLAGS += -DDEBUG -g 
endif 
else
CXX = cl.exe
CXX_FLAGS = /std:c++17 /DWIN32=1 /EHsc /Fe$@.exe
ifdef DEBUG 
CXX_FLAGS += /DDEBUG=1 /Zi
endif
endif

aoc:
	$(CXX) $(CXX_FLAGS) aoc.cc
