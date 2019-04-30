# Makefile for site_factor/xfactor plugin

AUTOMAKE_OPTIONS = foreign

PLUGIN_FLAGS = -module -avoid-version --export-dynamic

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src/common

pkglib_LTLIBRARIES = site_factor_xfactor.la

# site_factor NONE plugin.
site_factor_xfactor_la_SOURCES = site_factor_xfactor.c
site_factor_xfactor_la_LDFLAGS = $(PLUGIN_FLAGS) -lm

force:
$(site_factor_xfactor_la_LIBADD) : force
	@cd `dirname $@` && $(MAKE) `basename $@`
