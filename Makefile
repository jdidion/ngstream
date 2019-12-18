module = ngstream
repo = jdidion/$(module)
desc = Release $(version)
tests = tests
pytestops = -vv -s --full-trace

all: clean install test

build:
	python setup.py build_ext -i
	python setup.py sdist bdist_wheel

install: clean build
	python setup.py install $(installargs)

test:
	py.test $(pytestops) $(tests)

docs:
	make -C doc html

lint:
	pylint $(module)

clean:
	rm -Rf __pycache__
	rm -Rf **/__pycache__/*
	rm -Rf **/*.c
	rm -Rf **/*.so
	rm -Rf **/*.pyc
	rm -Rf dist
	rm -Rf build
	rm -Rf .adapters
	rm -Rf atropos.egg-info

docker:
	# build
	docker build -f Dockerfile -t $(repo):$(version) .
	# add alternate tags
	docker tag $(repo):$(version) $(repo):latest
	# push to Docker Hub
	docker login && docker push $(repo)

release:
	$(clean)
	# tag
	git tag $(version)
	# build
	$(BUILD)
	$(TEST)
	python setup.py sdist bdist_wheel
	# release
	python setup.py sdist upload -r pypi
	git push origin --tags
	$(github_release)
	$(docker)

tag:
	git tag $(version)

release: clean tag install test
	python setup.py sdist upload -r pypi
	git push origin --tags
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
