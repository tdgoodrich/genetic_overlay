clean_data:
	python Code/driver.py build_clean_data
arc:
	python Code/driver.py check_arcs
edgelist:
	python Code/driver.py build_edgelist
gephi:
	python Code/driver.py gephi
ppi:
	python Code/parse_ppi.py build_graph
correlation:
	python Code/driver.py correlation
