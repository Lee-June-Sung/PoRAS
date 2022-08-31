import sys, glob, os
import subprocess as sub

dir = sys.argv[1]



file_list = glob.glob('%s/*sam'%(dir))
file_list.sort()

w = open('Total_mapped_read.txt', 'w')

for file in file_list:
        print file
	cmd = 'samtools view -@ 10 -F 0x904 -c %s'%(file)
        print cmd
	p = sub.Popen(cmd.split(), stdout=sub.PIPE,stderr=sub.PIPE)
	output, errors = p.communicate()
	mapped_read = output.strip()

	w.write('%s\t%s\n'%(file, mapped_read))

w.close()






#print file_list





#samtools view -c -F 260 SAMPLE.bam

 #       p = sub.Popen(command_line7.split(), stdout=sub.PIPE,stderr=sub.PIPE)
  #      output, errors = p.communicate()
   #     mapped_read = output.strip()

