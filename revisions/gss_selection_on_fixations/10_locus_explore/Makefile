dbfiles= himu.1698870.sqlite3 himu.30884989.sqlite3 himu.33362904.sqlite3 himu.35868890.sqlite3 himu.37926691.sqlite3 himu.7234088.sqlite3 himu.8818318.sqlite3 himu.8968385.sqlite3 lowmu.16317964.sqlite3 lowmu.16631709.sqlite3 lowmu.18266900.sqlite3 lowmu.18378992.sqlite3 lowmu.22371142.sqlite3 lowmu.30377384.sqlite3 lowmu.35467076.sqlite3 lowmu.4215383.sqlite3 midmu.1231116.sqlite3 midmu.16185522.sqlite3 midmu.25492824.sqlite3 midmu.25938222.sqlite3 midmu.27933322.sqlite3 midmu.28769638.sqlite3 midmu.36197262.sqlite3 midmu.39937814.sqlite3

all: $(dbfiles)

%.sqlite3: %.gz process_population.py
	python3 process_population.py --infile $< --outfile $@
