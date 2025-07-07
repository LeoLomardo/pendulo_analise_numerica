all:
	gcc -o main main.c rk.c pendulo.c -lm
	./main > output/pendulo.csv
	python output/plot.py