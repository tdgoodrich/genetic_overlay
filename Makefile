clean_data:
	python Code/driver.py build_clean_data

arc:
	python Code/driver.py check_arcs

plot:
	python Code/driver.py plot_graph

gephi:
	python Code/driver.py gephi

edgelist:
	python Code/driver.py build_edgelist

ppi:
	python Code/parse_ppi.py build_graph

cluster:
	python Code/driver.py cluster
