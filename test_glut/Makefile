MAKE = mingw32-make --no-print-directory

SUBDIRS = cad_view msh_view field_evaluate scalar2d scalar3d solid2d solid3d fluid2d dkt2d helmholtz2d

all :
	$(MAKE) -C cad_view
	$(MAKE) -C msh_view
	$(MAKE) -C field_evaluate
	$(MAKE) -C scalar2d
	$(MAKE) -C scalar3d
	$(MAKE) -C solid2d
	$(MAKE) -C solid3d
	$(MAKE) -C fluid2d
	$(MAKE) -C dkt2d
	$(MAKE) -C helmholtz2d

clean :
	$(MAKE) clean -C cad_view
	$(MAKE) clean -C msh_view
	$(MAKE) clean -C field_evaluate
	$(MAKE) clean -C scalar2d
	$(MAKE) clean -C scalar3d
	$(MAKE) clean -C solid2d
	$(MAKE) clean -C solid3d
	$(MAKE) clean -C fluid2d
	$(MAKE) clean -C dkt2d
	$(MAKE) clean -C helmholtz2d


