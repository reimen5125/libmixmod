BUILD = build
BUILD_EX = EXAMPLES/build
SRC = SRC

INSTALL_PREFIX ?= $(HOME)/software

all: mixmod examples
.FORCE:

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD_EX):
	mkdir -p $(BUILD_EX)

mixmod: .FORCE | $(BUILD)
	cd $(BUILD) && cmake .. -DCMAKE_INSTALL_PREFIX=$(INSTALL_PREFIX)/libmixmod
	+$(MAKE) -C $(BUILD) install

examples: mixmod | $(BUILD_EX)
	cd $(BUILD_EX) && cmake .. -DCMAKE_INSTALL_PREFIX=$(INSTALL_PREFIX)/libmixmod
	+$(MAKE) -C $(BUILD_EX) install

clean:
	rm -rf $(BUILD) $(BUILD_EX)
