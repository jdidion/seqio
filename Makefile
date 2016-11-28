build = python setup.py install && \
		nose2 -C tests --coverage-report term-missing --coverage-config .coveragerc

install:
	$(call build,)

release:
	# tag
	git tag $(version)
	# build
	$(call build,)
	python setup.py sdist bdist_wheel
	# release
	twine register dist/seqio-$(version).tar.gz
	twine upload dist/seqio-$(version).tar.gz
	# push new tag after successful release
	git push origin --tags

docs:
	make -C docs api
	make -C docs html

readme:
	pandoc --from=markdown --to=rst --output=README.rst README.md
