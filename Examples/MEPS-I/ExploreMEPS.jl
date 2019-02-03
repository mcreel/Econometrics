# Read in the data

data = readdlm("meps1996.data")

#	The variables and their columns
#	VARIABLE						column
# 	Office based doctor visits	    1
# 	Outpatient doctor visits		2
# 	Emergency room visits			3
# 	Inpatient visits				4
# 	Dental visits					5
# 	Prescriptions					6
#  	PUBLIC_INS  					7
# 	PRIVATE_INS						8
# 	SEX								9
# 	AGE								10
# 	EDUC							11
# 	INCOME							12


names = ["OBDV", "OPV","IPV","ERV","DV","RX","PUB","PRIV","SEX","AGE","EDUC","INC/1000"]
data[:,12] = data[:,12]/1000.0
n = size(data,1)
println("MEPS data, 1996, complete data set statistics")
println("observations: ",n)
dstats(data, names);
""
