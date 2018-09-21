# Makefile for shared library

# C++ compiler
CXX = g++
# C++ flags
CXXFLAGS = -std=c++14 -fPIC -O3
# linking flags
LDFLAGS = -shared -ldcubic
# target lib
TARGET_LIB = lib/libdynamicboltz.so
# build dir
BUILD_DIR = build
# source dir
SOURCE_DIR = src
# install dir
INSTALL_LIB_DIR = /usr/local/lib
INSTALL_INCLUDE_DIR = /usr/local/include

# source files
SRC_NAMES = species.cpp \
	ixn_param.cpp \
	diff_eq_rhs.cpp \
	general.cpp \
	moment.cpp \
	unit.cpp \
	ixn_dicts.cpp \
	connections.cpp \
	lattice.cpp \
	adjoint.cpp
SRCS = $(addprefix $(SOURCE_DIR)/, $(SRC_NAMES))
OBJS = $(addprefix $(BUILD_DIR)/, $(SRC_NAMES:.cpp=.o))
DEPS = $(OBJS:.o=.d)

.PHONY: all clean install

all: $(BUILD_DIR) $(TARGET_LIB)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR) $(BUILD_DIR)/include_impl

$(TARGET_LIB): $(OBJS)
	$(CXX) ${LDFLAGS} $^ -o $@

# Every .d file will be matched here, and will be compiled from a .cpp file
$(OBJS): $(BUILD_DIR)/%.o : $(SOURCE_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

# Use of '-' in front will ignore the errors since the files arent found right away
ifneq ($(MAKECMDGOALS),clean)
-include $(DEPS)
endif

clean:
	$(RM) ${TARGET_LIB} $(OBJS) $(DEPS)

install:
	mkdir -p $(INSTALL_LIB_DIR)
	mkdir -p $(INSTALL_INCLUDE_DIR)
	mkdir -p $(INSTALL_INCLUDE_DIR)/dynamicboltz_bits
	mkdir -p $(INSTALL_INCLUDE_DIR)/dynamicboltz_bits/fwds
	cp -p $(TARGET_LIB) $(INSTALL_LIB_DIR)
	cp -p include/dynamicboltz_bits/*.hpp $(INSTALL_INCLUDE_DIR)/dynamicboltz_bits
	cp -p include/dynamicboltz_bits/fwds/*.hpp $(INSTALL_INCLUDE_DIR)/dynamicboltz_bits/fwds
	cp -p include/dynamicboltz $(INSTALL_INCLUDE_DIR)
