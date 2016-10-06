
ifndef LARLITECV_BASEDIR
ERROR_MESSAGE := $(error LARLITECV_BASEDIR is not set... run configure.sh!)
endif

ifndef LARLITE_BASEDIR
ERROR_MESSAGE := $(error LARLITE_BASEDIR is not set... run configure.sh!)
endif

ifndef LARCV_BASEDIR
ERROR_MESSAGE := $(error LARCV_BASEDIR is not set... run configure.sh!)
endif

OSNAME          = $(shell uname -s)
HOST            = $(shell uname -n)
OSNAMEMODE      = $(OSNAME)

include $(LARLITECV_BASEDIR)/Makefile/Makefile.${OSNAME}

CORE_SUBDIRS := Hashlib2plus Base
#ifeq ($(LARLITECV_NUMPY),1)
#CORE_SUBDIRS += PyUtil
#endif
#ifeq ($(LARLITECV_OPENCV),1)
#  CORE_SUBDIRS += CVUtil
#endif

APP_SUBDIRS := ThruMu

.phony: all clean

all: obj lib

clean: clean_app clean_core
	@rm -f $(LARLITECV_LIBDIR)/liblarlitecv.so
clean_core:
	@for i in $(CORE_SUBDIRS); do ( echo "" && echo "Cleaning $$i..." && cd $(LARLITECV_COREDIR)/$$i && rm -rf $(LARLITECV_BUILDDIR)/$$i && rm -rf $(LARLITECV_BUILDDIR)/lib/*$ii.* ) || exit $$?; done
clean_app:
	@for i in $(APP_SUBDIRS); do ( echo "" && echo "Cleaning $$i..." && cd $(LARLITECV_APPDIR)/$$i && rm -rf $(LARLITECV_BUILDDIR)/$$i && rm -rf $(LARLITECV_BUILDDIR)/lib/*$ii.* ) || exit $$?; done

obj:
	@echo
	@echo Building core...
	@echo
	@for i in $(CORE_SUBDIRS); do ( echo "" && echo "Compiling $$i..." && cd $(LARLITECV_COREDIR)/$$i && $(MAKE) ) || exit $$?; done
	@echo Building app...
	@for i in $(APP_SUBDIRS); do ( echo "" && echo "Compiling $$i..." && cd $(LARLITECV_APPDIR)/$$i && $(MAKE) ) || exit $$?; done

lib: obj
	@ echo
	@ if [ `python ${LARLITECV_BASEDIR}/bin/libarg.py build` ]; then \
	    echo Linking library...; \
	    $(SOMAKER) $(SOFLAGS) $(shell python $(LARLITECV_BASEDIR)/bin/libarg.py); \
	  else \
	   echo Nothing to be done for lib...; \
	fi
	@echo 
