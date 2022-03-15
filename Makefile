ALL:
	number=1 ; while [[ $$number -le 15 ]] ; do \
		echo "Number: $$number" ; \
		num1=`echo $${number}|awk '{printf("%02d\n",$$0)}'`; \
		echo "num1: $$num1" ; \
		xelatex -shell-escape  $${num1}-discussion.tex; \
		bibtex  $${num1}-discussion; \
		xelatex  $${num1}-discussion.tex; \
		((number = number + 1)) ; \
	done
	