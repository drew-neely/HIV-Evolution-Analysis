import csv
from os import listdir
from os.path import isfile, join
from itertools import combinations
from re import sub

data_index_filename = 'data_index.csv'
ggdc_data_directory = 'ggdc_data'
data_column_names = ["query", "reference", "f1-ddh", "f1-ci", "f1-dist", "f1-prob", 
						"f2-ddh", "f2-ci", "f2-dist", "f2-prob", 
						"f3-ddh", "f3-ci", "f3-dist", "f3-prob", "gc-diff"]
to_filename = lambda str : "ggdc_results_" + str + ".csv"

class DisimilarityMatrix :

	def __init__(self, sample_ids) :
		self.sample_ids = sample_ids
		self.matrix = {pair: None for pair in combinations(self.sample_ids, 2)}
		
	def add(self, sample_pair, dist) :
		key = sample_pair if sample_pair in self.matrix else (sample_pair[1], sample_pair[0])
		assert key in self.matrix, "Invalid pair passed to DisimilarityMatrix.add"
		assert not self.matrix[key], "Attempting to overwrite disimilarity matrix entry"
		self.matrix[key] = dist
	
	def __str__(self) :
		res = ""
		for y in range(0, len(self.sample_ids)) :
			for x in range(0, y) :
				pair = (self.sample_ids[y], self.sample_ids[x])
				pair = pair if pair in self.matrix else (pair[1], pair[0])
				assert pair in self.matrix, "Invalid pair in DisimilarityMatrix.__str__"
				res += self.matrix[pair] + " "
			res += "0\n"
		return res


# Read in data_index_filename
# fields: 
	# 	patient_id
	# 	days_from_seroconversion
	# 	sample_id (This is "<patient_id>-<days_from_seroconversion>")
	# 	accession_ids
	# 	job_id
file = open(data_index_filename, mode='r')
index = csv.DictReader(file)
index = [r for r in index] # Dict readaer is dumb and goes away once you read from it
file.close()
job_file_names = set([to_filename(r["job_id"]) for r in index if r["job_id"] != "NA"])

# Get ggdc_data file name list
ggdc_data_filenames = set([f for f in listdir(ggdc_data_directory) if isfile(join(ggdc_data_directory, f))])
ggdc_data_filenames = ggdc_data_filenames.intersection(job_file_names)
if len(ggdc_data_filenames) != len(job_file_names) :
	missing = list(job_file_names.difference(ggdc_data_filenames))
	print("Missing the following expected files:")
	print(missing)
ggdc_data_filenames = [join(ggdc_data_directory, name) for name in ggdc_data_filenames]

# init similarity matrix
mat = DisimilarityMatrix([r["sample_id"] for r in index])

# init lookup table (accesion -> sample_id)
sample_ids = {}
for row in index:
	key = row["accession_ids"].strip().replace("\s+", " ")
	sample_ids[key] = row["sample_id"]

# iterate over data files
for data_filename in ggdc_data_filenames :
	with open(data_filename, mode='r') as file:
		data = csv.DictReader(file, fieldnames=data_column_names)
		# ignore first two rows on account of these csv files being dumb
		data = list(data)[2:]
		for d in data :
			qry_asc = sub("\s+", " ", d['query'].strip())
			ref_asc = sub("\s+", " ", d['reference'].strip())
			assert qry_asc in sample_ids, "invalid query: " + qry_asc + " in " + data_filename
			assert ref_asc in sample_ids, "invalid reference: " + ref_asc + " in " + data_filename
			q_id = sample_ids[qry_asc]
			r_id = sample_ids[ref_asc]
			mat.add((q_id, r_id), d["f2-dist"])

with open("disimilarity_matrix.txt", "w") as outfile :
	outfile.write(str(mat))
			


	