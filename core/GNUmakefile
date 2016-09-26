
ifndef LARLITECV_COREDIR
ERROR_MESSAGE := $(error LARLITECV_COREDIR is not set... run configure.sh!)
endif

#####################################################################################
#
# Define directories to be compile upon a global "make"...
#
SUBDIRS := Hashlib2plus Base

#####################################################################################
#
# COMPILATION...
#
#.phony: all configure default-config clean
.phony: all clean

all:
	@for i in $(SUBDIRS); do ( echo "" && echo "Compiling $$i..." && cd $$i && $(MAKE) ) || exit $$?; done
#####################################################################################
#
# CLEANs...
#
clean:
	@for i in $(SUBDIRS); do ( echo "" && echo "Cleaning $$i..." && cd $$i && $(MAKE) clean && rm -rf $(LARLITECV_BUILDDIR)/$$i && rm -rf $(LARLITECV_BUILDDIR)/lib/*$ii.* ) || exit $$?; done

#####################################################################################
#
# DOCUMENTATION...
#
doxygen:
	@echo 'dOxygenising LArLiteCV core...'
	@mkdir -p $(LARLITECV_COREDIR)/doc/dOxygenMyProject
	@doxygen $(LARLITECV_COREDIR)/doc/doxygenMyProject.script
#
#####################################################################################
#EOF
