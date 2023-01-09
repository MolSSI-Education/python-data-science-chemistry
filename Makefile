
.PHONY: book 

html:
	jupyter-book build book

clean: book/_build
	echo "Removing everything under _build"
	rm -rf book/_build
