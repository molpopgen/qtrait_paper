all: pop.gz pop.sqlite3 pop2.gz pop2.sqlite3 pop3.gz pop3.sqlite3 StrengthOfSelection.pdf

pop.gz: simulate.py
	python3 simulate.py --popsize 5000 --filename pop.gz

pop.sqlite3: process_population.py pop.gz
	python3 process_population.py --infile pop.gz --outfile pop.sqlite3

pop2.gz: simulate.py
	python3 simulate.py --popsize 5000 --filename pop2.gz --seed 11018

pop2.sqlite3: process_population.py pop2.gz
	python3 process_population.py --infile pop2.gz --outfile pop2.sqlite3

pop3.gz: simulate.py
	python3 simulate.py --popsize 5000 --filename pop3.gz --seed 27371

pop3.sqlite3: process_population.py pop3.gz
	python3 process_population.py --infile pop3.gz --outfile pop3.sqlite3

StrengthOfSelection.pdf: pop.sqlite3 pop2.sqlite3 pop3.sqlite3 plot.py
	python3 plot.py
