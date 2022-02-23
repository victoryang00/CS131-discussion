ALL:
	xelatex -shell-escape  02-discussion.tex
	bibtex  02-discussion
	xelatex  02-discussion.tex