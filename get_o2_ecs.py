from Bio.ExPASy import Enzyme

infile = "/Users/jagodajablonska/oxyphen/DATA/enzyme.dat"
handle = open(infile)
records = Enzyme.parse(handle)

O2_ec_list = []
existing_list  =  open("/Users/jagodajablonska/oxyphen/DATA/oxygen_ecclasses").read().splitlines()


for record in records:
	EC_num = record['ID']
	reaction = record['CA']

	if "=" in reaction:
		substrates = [x.strip() for x in reaction.split("=")[0].split("+")]
		products = [x.strip() for x in reaction.split("=")[1].split("+")]

	if "O(2)" in substrates or "O(2)" in products:
		if EC_num not in existing_list:
			print(EC_num)



