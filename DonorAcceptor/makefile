# Compiler
CXX = mpicxx

# Compiler flags
CXXFLAGS = -Wall -O2 -I.

# Include directories (if needed)
INCLUDES = -I.

# Source files
SRCS = DataStruct.cpp IO.cpp ReadGrid.cpp WriteGrid.cpp FindFringePoints.cpp main.cpp

# Object directory
OBJDIR = obj

# Object files (stored in obj directory)
OBJS = $(patsubst %.cpp,$(OBJDIR)/%.o,$(SRCS))

# Executable name
EXEC = overset

# Default target
all: $(EXEC)

# Create the object directory
$(OBJDIR):
	mkdir -p $(OBJDIR)

# Linking the executable
$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Compile each .cpp file into .o and store in obj directory
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean up build files
clean:
	rm -rf $(OBJDIR) $(EXEC)
