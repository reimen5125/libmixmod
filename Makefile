BUILD = build
BUILD_EX = EXAMPLES/build
SRC = SRC

DEBUG ?= 0
INSTALL_PREFIX ?= $(HOME)/software

CMAKEFLAGS = -DCMAKE_INSTALL_PREFIX=$(INSTALL_PREFIX)/libmixmod

ifeq ($(DEBUG), 1)
CMAKEFLAGS += -DCMAKE_BUILD_TYPE="Debug"
else
CMAKEFLAGS += -DCMAKE_BUILD_TYPE="Release"
endif

all: mixmod examples
.FORCE:

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD_EX):
	mkdir -p $(BUILD_EX)

mixmod: .FORCE | $(BUILD)
	cd $(BUILD) && cmake .. $(CMAKEFLAGS)
	+$(MAKE) -C $(BUILD) install

examples: mixmod | $(BUILD_EX)
	cd $(BUILD_EX) && cmake .. $(CMAKEFLAGS)
	+$(MAKE) -C $(BUILD_EX) install

clean:
	rm -rf $(BUILD) $(BUILD_EX)
