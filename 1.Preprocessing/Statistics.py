import sys, glob, os
#from tqdm import tqdm

#in_dir = sys.argv[1]
in_dir = os.getcwd() + '/'

file_list = glob.glob('%s/*.f*q'%(in_dir))
file_list.sort()

statistics = 'Sample' +'\t'+ 'Reads' +'\t'+ 'Total_length' +'\t'+ 'read_average' +'\t'+ 'Q20' +'\t'+ 'Q30'

count = 0

for name in file_list:
	file = open(name, 'r')
	a = 0
	Reads = 0
	Total_length = 0
	Q20_count = 0
	Q30_count = 0

	count += 1
	print '%s/%s Read_statistics running...'%(count, len(file_list))

	for line in file.readlines():
		if line == '+\n':
			a = 1
			continue

		if line[0:4] == '+SRR':
                        a = 1
                        continue

                if line[0:4] == '+ERR':
                        a = 1
                        continue

                if line[0:4] == '+DRR':
                        a = 1
                        continue


		if a == 1:
			a = 0
			Reads += 1
			Total_length += len(line.strip())
			Q20_count += line.count('5') + line.count('6') + line.count('7') + line.count('8') +line.count('9') + line.count(':') + line.count(';') + line.count('<') + line.count('=') + line.count('>') + line.count('?') + line.count('@') + line.count('A') + line.count('B') + line.count('C') + line.count('D') + line.count('E') + line.count('F') + line.count('G') + line.count('H') + line.count('I') + line.count('J') + line.count('K')
			Q30_count += line.count('?') + line.count('@') + line.count('A') + line.count('B') + line.count('C') + line.count('D') + line.count('E') + line.count('F') + line.count('G') + line.count('H') + line.count('I') + line.count('J') + line.count('K')

	read_average = float(Total_length)/float(Reads)
	Q20 = float(Q20_count)/float(Total_length)*100
	Q30 = float(Q30_count)/float(Total_length)*100
	Sample = name.strip().split('/')[-1]

	statistics += '\n' + Sample +'\t'+ str(Reads) +'\t'+ str(Total_length) +'\t'+ str(read_average) +'\t'+ str(Q20) +'\t'+ str(Q30)

file1 = open('%sRead_statistics'%(in_dir), 'w')
file1.write(statistics)
file1.close()

print '\n%s/%s Read_statistics Finish!!'%(len(file_list), len(file_list))
