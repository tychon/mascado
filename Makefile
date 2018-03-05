
test:
	python -m unittest discover --start-directory tests tests

docupdate:
	rm -r docs/source/maskastrometry.*
	cd docs; \
	sphinx-apidoc -MeT -o source/ ../maskastrometry

docbuild:
	cd docs; \
	make html

docshow:
	xdg-open docs/build/html/py-modindex.html
