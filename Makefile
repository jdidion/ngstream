package = ngstream
repo = jdidion/$(package)
desc = Release $(version)
tests = tests
pytestops = -vv -s --full-trace

all: clean install install_extra_requirements install_test_requirements test test_release_setup

build:
	python setup.py build_ext -i
	python setup.py sdist bdist_wheel

install: clean build
	pip install --upgrade dist/*.whl $(installargs)

install_test_requirements:
	pip install -r requirements-test.txt

install_extra_requirements:
	pip install -r requirements-extra.txt

test: install install_extra_requirements install_test_requirements
	pytest $(pytestops) $(tests)

test_release_setup:
	twine check dist/*

docs:
	make -C doc html

lint:
	pylint $(package)

reformat:
	black $(package)
	black $(tests)

clean:
	rm -Rf __pycache__
	rm -Rf **/__pycache__/*
	rm -Rf **/*.c
	rm -Rf **/*.so
	rm -Rf **/*.pyc
	rm -Rf dist
	rm -Rf build
	rm -Rf .adapters
	rm -Rf $(module).egg-info

tag:
	git tag $(version)


push_tag:
	git push origin --tags

del_tag:
	git tag -d $(version)

pypi_release:
	twine upload dist/*

release: clean tag
	${MAKE} test pypi_release push_tag || (${MAKE} del_tag && exit 1)

	# github release
	curl -v -i -X POST \
		-H "Content-Type:application/json" \
		-H "Authorization: token $(token)" \
		https://api.github.com/repos/$(repo)/releases \
		-d '{\
		  "tag_name":"$(version)",\
		  "target_commitish": "master",\
		  "name": "$(version)",\
		  "body": "$(desc)",\
		  "draft": false,\
		  "prerelease": false \
		}'
