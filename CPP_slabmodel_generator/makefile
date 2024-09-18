# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Werror -Wextra -std=c++17 -O2

# Directories
SRC_DIRS = . Molecule HITRAN
OBJ_DIR = obj

# Source files
SRCS = $(foreach dir, $(SRC_DIRS), $(wildcard $(dir)/*.cc))

# Object files
OBJS = $(patsubst %.cc, $(OBJ_DIR)/%.o, $(SRCS))

# Executable name
TARGET = main

# Default target to compile app
all: $(TARGET)

# Linking object files to create executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

# Rule for compiling .cc files into .o files
$(OBJ_DIR)/%.o: %.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o:: %/%.cc
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean