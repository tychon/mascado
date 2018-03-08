
test:
	python -m unittest discover --start-directory tests

docupdate:
	rm -rf docs/source/mascado.*
	cd docs; \
	sphinx-apidoc -MeT -o source/ ../mascado

docbuild:
	cd docs; \
	make html

docshow:
	xdg-open docs/build/html/py-modindex.html
