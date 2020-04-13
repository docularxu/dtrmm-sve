ifeq ($(BLISLAB_USE_INTEL),true)
include $(BLISLAB_DIR)/make.inc.files/make.intel.inc
else
include $(BLISLAB_DIR)/make.inc.files/make.gnu.inc
endif

$(info * Using CFLAGS=${CFLAGS})
$(info * Using LDFLAGS=${LDFLAGS})

FRAME_CC_SRC=  \
							 dtrmm/my_dtrmm.c \
							 dtrmm/bl_dtrmm_ref.c \
							 dtrmm/bl_dgemm_util.c \

FRAME_CPP_SRC= \

KERNEL_SRC=    \
							 kernels/bl_dtrmm_ukr.c \

OTHER_DEP = \
			                 include/bl_dtrmm.h \
			                 include/bl_dgemm_kernel.h \

BLISLAB_OBJ=$(FRAME_CC_SRC:.c=.o) $(FRAME_CPP_SRC:.cpp=.o) $(KERNEL_SRC:.c=.o) $(FRAME_CC_SRC_S:.c=.os) $(KERNEL_SRC_S:.c=.os)

all: $(LIBBLISLAB) $(SHAREDLIBBLISLAB) TESTBLISLAB

TESTBLISLAB: $(LIBBLISLAB)
	cd $(BLISLAB_DIR)/test && $(MAKE) && cd $(BLISLAB_DIR)

$(LIBBLISLAB): $(BLISLAB_OBJ)
	$(ARCH) $(ARCHFLAGS) $@ $(BLISLAB_OBJ)
	$(RANLIB) $@

$(SHAREDLIBBLISLAB): $(BLISLAB_OBJ)
	$(CC) $(CFLAGS) -shared -o $@ $(BLISLAB_OBJ) $(LDLIBS)

# ---------------------------------------------------------------------------
# Object files compiling rules
# ---------------------------------------------------------------------------
%.o: %.c $(OTHER_DEP)
	$(CC) $(CFLAGS) -c $< -o $@ $(LDFLAGS)

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $< -o $@ $(LDFLAGS)
# ---------------------------------------------------------------------------

clean:
	-rm $(BLISLAB_OBJ) $(LIBBLISLAB) $(SHAREDLIBBLISLAB) dtrmm/*~ kernels/*~ kernels/*.o test/*~ include/*~ *~ make.inc.files/*~
	$(MAKE) clean -C test

