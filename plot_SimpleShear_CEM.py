#!/usr/bin/python

def trova_punto (nome_file_sorgente, elenco_picchi_CEM): # nome_file_sorgente da elenco_sorgenti

	file_outcomes = open (nome_file_sorgente, 'r')
	file_outcomes_lines = []
	for line in file_outcomes.readlines():
		file_outcomes_lines.append(line.rstrip('\n'))
	file_outcomes.close()

	file_outcomes_dictionary = {} # {genome_location_string : n_of_distinct_shear_float}
	total_n_of_distinct_shear = 0
	for line in file_outcomes_lines[1:]:
		line_split = line.split('\t')
		key = line_split[0] + " " + line_split[1] + " " + line_split[2]
		n_of_distinct_shear = int(line_split[7])
		total_n_of_distinct_shear = total_n_of_distinct_shear + n_of_distinct_shear # Normalization factor
		file_outcomes_dictionary.update({key:n_of_distinct_shear})

	cem_dictionary = {}
	for cem in elenco_picchi_CEM:
		cem_split = cem.split('\t')
		genome_location = cem_split[0]
		#cem_symbol = cem_split[1]
		n_of_distinct_shear = 0
		if (file_outcomes_dictionary.has_key(genome_location)):
			n_of_distinct_shear = n_of_distinct_shear + file_outcomes_dictionary[genome_location]
		#account for "shoulders"
		genome_location_split = genome_location.split(' ')
		genome_location_back = genome_location_split[0] + " " + str(int(genome_location_split[1]) - 1) + " " + genome_location_split[2]
		print genome_location_back
		genome_location_forward = genome_location_split[0] + " " + str(int(genome_location_split[1]) + 1) + " " + genome_location_split[2]
		print genome_location_forward
		if (file_outcomes_dictionary.has_key(genome_location_back)):
			n_of_distinct_shear = n_of_distinct_shear + file_outcomes_dictionary[genome_location_back]
		if (file_outcomes_dictionary.has_key(genome_location_forward)):
			n_of_distinct_shear = n_of_distinct_shear + file_outcomes_dictionary[genome_location_forward]
		n_of_distinct_shear = float(n_of_distinct_shear) / float(total_n_of_distinct_shear) # Normalize n_of_distinct_shear
		cem_dictionary.update({genome_location:n_of_distinct_shear})

	nome_file_sorgente_split = nome_file_sorgente.split('.')
	diluizione = nome_file_sorgente_split[0] + "\n"

	###Changes to plot###
	if ((diluizione == "A\n") or (diluizione == "M\n")):
		diluizione = "100.0"
	elif ((diluizione == "B\n") or (diluizione == "N\n")):
		diluizione = "75.0"
	elif ((diluizione == "C\n") or (diluizione == "O\n")):
		diluizione = "50.0"
	elif ((diluizione == "D\n") or (diluizione == "P\n")):
		diluizione = "25.0"
	elif ((diluizione == "E\n") or (diluizione == "Q\n")):
		diluizione = "10.0"
	elif ((diluizione == "F\n") or (diluizione == "R\n")):
		diluizione = "1.0"
	elif ((diluizione == "G\n") or (diluizione == "S\n")):
		diluizione = "0.1"
	elif ((diluizione == "H\n") or (diluizione == "T\n")):
		diluizione = "0.01"
	elif ((diluizione == "I\n") or (diluizione == "U\n")):
		diluizione = "0.001"
	elif ((diluizione == "L\n") or (diluizione == "V\n")):
		diluizione = "0.000"
	#####################

	lista_abundance_cem = []
	ordered_genome_locations = sorted(cem_dictionary.keys())
	for key in ordered_genome_locations:
		lista_abundance_cem.append("\t" + str(cem_dictionary[key]))
	lista_abundance_cem.append("\n")

	return diluizione, lista_abundance_cem



# Elenco dei file sorgenti "*outcomes.tsv"
file_elenco_sorgenti = open ('file_sorgenti_gatc_FB386388_p20_f50ec1o3.txt', 'r') #Anche con file_sorgenti_gsk_FB386388_p20_f50ec1o3.txt
#file_elenco_sorgenti = open ('file_sorgenti_gsk_FB386388_p20_f50ec1o3.txt', 'r')
elenco_sorgenti = []
for line in file_elenco_sorgenti.readlines():
	elenco_sorgenti.append(line.rstrip('\n'))
file_elenco_sorgenti.close()

# Elenco dei picchi CEM
file_CEM = open ('sorgente_CEM.txt')
elenco_picchi_CEM = []
for line in file_CEM.readlines():
	elenco_picchi_CEM.append(line.rstrip('\n'))
file_CEM.close()

#File output
file_output = open ('gatc_plot_SimpleShear.txt', 'w') #Da cambiare
#file_output = open ('gsk_plot_SimpleShear.txt', 'w')

for nome_file_sorgente in elenco_sorgenti:
	diluizione, lista_abundance_cem = trova_punto (nome_file_sorgente, elenco_picchi_CEM)
	file_output.write(diluizione)
	file_output.writelines(lista_abundance_cem)
file_output.close()