SHELL := /bin/bash
ALL:
	number=1 ; while [[ $$number -le 12 ]] ; do \
		echo "Number: $$number" ; \
		num1=`echo $${number}|awk '{printf("%02d\n",$$0)}'`; \
		echo "num1: $$num1" ; \
		xelatex -interaction=nonstopmode -shell-escape $${num1}-discussion.tex; \
		bibtex  $${num1}-discussion; \
		xelatex -interaction=nonstopmode -shell-escape $${num1}-discussion.tex; \
		((number = number + 1)) ; \
	done
	
