#!C:\Python26\python.exe
#useful functions for bioinformatics work
#jqian
#Aug 5th,2010
#kgi/TMP room B213
#all tested if no special notes
#################################
#GC_content(dna) this is to calculate the gc content of a dna sequence
#complement(dna) this is to do the complementary sequence for dna
#restrict(dna,  enz) Restriction site occurrences as a list (),find  all  start  positions  of  a  restriction  site
#codons(s,frame=0) Get the codon list from a DNA sequence ()
#revcomp(dna) #Reverse Complement of DNA,reverse  complement  of  a  DNA  sequence
#################################
import todo
import elixir
import configure_finger
import TmDeltaG
import SecStructures_jf4
from Bio.SeqRecord import SeqRecord
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
#################################
#this is to calculate the gc content of a dna sequence
import os
import SeqDep
from string import *
def GC_content(dna):
    gc  =  (count(dna,  'C')  +  count(dna,  'c')  +  count(dna,  'G')  +  count(dna,  'g'))  /  float(len(dna))  *  100
    return gc

###########################################
#this is to do the complementary sequence for dna
import string
def complement(dna):
    dna = dna.encode('ascii')
    t=string.maketrans("AGCTagct",	"TCGAtcga") # to make the translation table
    return translate(dna, t)

################################
#Restriction site occurrences as a list (),find  all  start  positions  of  a  restriction  site
from string import *
def  restrict(dna,  enz):
    res  =  []
    site = dna.find(enz)
    while (site != -1):
        res.append(site)
        site  =  dna.find(enz,  site  +  1)
    return  res

#####################################
#Get the codon list from a DNA sequence ()
from string import *
def codons(s,frame=0):
#frame can be 0, 1 or 2
    codons=[]
    end=len(s[frame:]) - (len(s[frame:]) % 3) - 1
    for i in range(frame,end,3):
        codons.append(s[i:i+3])
    return codons

##########################################
#Reverse Complement of DNA,reverse  complement  of  a  DNA  sequence
import  string
def  revcomp(dna):
    dna = dna.encode('ascii')
    comp = dna.translate(string.maketrans("AGCTagct", "TCGAtcga"))
    lcomp = list(comp)
    lcomp.reverse()
    return  join(lcomp,  "")

#################################
#write
#not tested yet
#from Bio import SeqIO
#input_file = open('NC_000962.fasta', 'r')
#output_file = open('NC_list.txt','wb')
#for cur_record in SeqIO.parse(input_file, "fasta") :
#         #output_file.write(gene_seq)
#         #print >>output_file, cur_record.seq
#         list1=restrict(cur_record.seq,  'GAGTC')
#         list1.sort()
#         list2=restrict(cur_record.seq,  'GACTC')
#         list2.sort()
#
#
#output_file.close()
#input_file.close()
############################
def site_find_finger(list1, list2):
#to find the location for all the recognition site
    pair=[]
    list1.sort()
    list2.sort()
    for i in range(len(list1)-1):
        if ((list1[i+1]-list1[i])<=39 and (list1[i+1]-list1[i])>=17): # this is HTT
            pair.append(list1[i])
            pair.append(list1[i+1])

    for i in range(len(list2)-1):
        if ((list2[i+1]-list2[i])<=39 and (list2[i+1]-list2[i])>=17): # this is TTH
            pair.append(list2[i])
            pair.append(list2[i+1])

    for i in range(len(list1)):
        for j in range(len(list2)):
            if (list2[j]-list1[i]<=43) and (list2[j]-list1[i]>=21) and (list2[j]<list1[i+1] or i==len(list2)) and (list2[j-1]<list1[i] or j==0): #this is HTH
                pair.append(list1[i])
                pair.append(list2[j])



    return pair

##########################
def find_finger_seq(seq,list1, list2,file, maximum=30,minimum=8,HTH=1,HTT=1,TTH=1):
#to find the location for all the recognition site and write the trigger to file
	list1.sort()
	list2.sort()
	rank=1
	outfile=open(file,'w')
	if (HTT==1):
	    for i in range(len(list1)-1):
		if ((list1[i+1]-list1[i])<=(int(maximum,10)+9) and (list1[i+1]-list1[i])>=(int(minimum,10)+9)):# this is HTT
		    location=list1[i]+1
		    length=list1[i+1]-list1[i]+9
		    finger=seq[list1[i]:list1[i+1]+9]
		    finger_inside=seq[list1[i]+9:list1[i+1]]
		    if finger_inside.find('GAGTC')==-1 and finger_inside.find('GACTC')==-1:
			print >>outfile,">seq_%s_HTT_%s_%s" %(rank,location,length)
			print >>outfile,finger
			rank=rank+1


	if (TTH==1):
	    for i in range(len(list2)-1):
		if ((list2[i+1]-list2[i])<=int(maximum,10)+9 and (list2[i+1]-list2[i])>=int(minimum,10)+9): # this is TTH
                        location=list2[i]+1-4
			length=list2[i+1]-list2[i]+9
			finger=seq[list2[i]-4:list2[i+1]+5]
			finger_inside=seq[list2[i]+5:list2[i+1]-4]
			if finger_inside.find('GAGTC')==-1 and finger_inside.find('GACTC')==-1:
				print >>outfile,">seq_%s_TTH_%s_%s" %(rank,location,length)
				print >>outfile,finger
				rank=rank+1


	if (HTH==1):
	    p=0
	    for i in range(len(list1)):
		for j in range(p,len(list2)):
                    if list2[j]<list1[i]:
			p=j
                        continue

                    elif list2[j]-list1[i]>int(maximum,10)+13 :
                        break

                    elif (list2[j]-list1[i]<=int(maximum,10)+13) and (list2[j]-list1[i]>=int(minimum,10)+13) and (i==len(list1)-1 or list2[j]<list1[i+1]) and (j==0 or list2[j-1]<list1[i]): #this is HTH
                                location=list1[i]+1
				length=list2[j]-list1[i]+5
				finger_inside=seq[list1[i]+5:list2[j]]
				if finger_inside.find('GAGTC')==-1 and finger_inside.find('GACTC')==-1:
				    print >>outfile,">seq_%s_HTH_%s_%s" %(rank,location,length)
				    print >>outfile,seq[list1[i]:list2[j]+5]
				    rank=rank+1


	outfile.close()

########################


########################
#get a sequence by accname from NCBI and write in in fasta format then store it in a file

def get_seq(accname,file):
    Entrez.email = 'jqian@kgi.edu'
    handle = Entrez.efetch(db='nucleotide',id=accname, rettype='fasta')
    record = SeqIO.read(handle, 'fasta')
    output_file = open(file,'wb')
    print >>output_file, record.format("fasta")
    handle.close()
    output_file.close()


#####################################
class Finger():
	def __init__(self,text):
		text = text.replace("\r\n","\n") #Crude way of dealing with \r\n
		assert text[0] == ">", text
		text = text.split("\n>",1)[0] # Only do the first record if more than one
		title, sequence = text.split("\n", 1)
		title = title[1:]
		self.title = title
		self.sequence = sequence.replace("\n","")
		temp= title.split("_")
		self.rank=temp[1]
		self.type=temp[2]
		self.len=len(self.sequence)
		#self.location=

########################################
##############################################
#read in the finger's name
def read_name(fasta_file):
    from Bio import SeqIO
    typ=[]
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        typ.append(seq_record.id)
    return typ


##############################################
#read in a dictionary between sequence name and rank
def read_diction_name(fasta_file):
    from Bio import SeqIO
    typ=[]
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
	temp=split(seq_record.id,"_")
	rank=temp[1]
	typ[seq_record.id]=rank

    return typ

##############################################
#read a list file
def read_list(listfile):
	f = open(listfile,'r')
	typ=[]
	while True:
		line = f.readline()
		line = line.rstrip('\n')
		if len(line) == 0: # Zero length indicates EOF
			break
	typ.append(line) #read one item for the list
	return typ
	f.close()

##############################################
#this is to see if the sequence is conserved in other species, the seq_name_list has some special requiremennt
from Bio.Blast import NCBIXML
from numpy import *
from string import *
def conservation_check(blast_result_file,conservation_list,seq_name_list):
     results = [[0] * len(conservation_list) for row in range(len(seq_name_list))]
     for record in NCBIXML.parse(open(blast_result_file)) :
          if record.alignments :
               name=record.query
               num=split(record.query,"_")[1]
               type=split(record.query,"_")[2]
               num=int(num,10)
               for align in record.alignments :
                    for hsp in align.hsps :
                         for index in range(len(conservation_list)):
                              if align.hit_def.find(conservation_list[index]) >-1:
                                   s1=bool(hsp.query_start ==1)
                                   s2=bool(hsp.query_start <= record.query_length-21)
                                   e1=bool(hsp.query_end == record.query_length)
                                   e2=bool(hsp.query_end >= 21)
                                   #print "result",results[num-1][index]
                                   #print "result",len(conservation_list)
                                   if (type =='HTH'):
                                        if (s2 and e1) and (results[num-1][index]==-1):
                                             results[num-1][index]=2
                                                #print "itshouldbe","num",num,"index",index,"change?",results[num][index+1]
                                        elif (s1  and e2) and (results[num-1][index]==1):
                                                results[num-1][index]=2
                                        elif (s1  and e2) and (results[num-1][index]==0):
                                                results[num-1][index]=-1
                                        elif (s2 and e1) and (results[num-1][index]==0):
                                                results[num-1][index]=1
						#print "into",num,index,results[num][index+1]
					    #print "HTH"
					    #print "hsp.query_start",hsp.query_start, "end",hsp.query_end, "len",record.query_length, "num",num, "index",index
					    #print "firt", bool(s2)
					    #print "second", bool(e1)
                                   elif (type =='HTT'):
					    #print "HTT",bool(s2 and e1)
                                        if(s2 and e1):
                                             results[num-1][index]=1
                                             #print "into",num,index,results[num-1][index]
					    #print "HTT"
                                   elif (type =='TTH'):
                                        if(s1  and e2):
						#print "into only ",num,index,name
                                             results[num-1][index]=-1
					    #print "TTH"
     print "doneconservation"
     return results




##########################################################
#this is to see if the sequence is conserved in other species, the seq_name_list has some special requiremennt
from Bio.Blast import NCBIXML
from numpy import *
from string import *
def conservation_check2(blast_result_file,conservation_list,seq_name_list):
     results = [[0] * len(conservation_list) for row in range(len(seq_name_list))]
     for record in NCBIXML.parse(open(blast_result_file)) :
          if record.alignments :
               name=record.query
               num=split(record.query,"_")[1]
               type=split(record.query,"_")[2]
               num=int(num,10)
               for align in record.alignments :
                    for hsp in align.hsps :
                         for index in range(len(conservation_list)):
                              if align.hit_def.find(conservation_list[index]) >-1:
                                   s1=bool(hsp.query_start ==1)
                                   s2=bool(hsp.query_start <=5)
                                   e1=bool(hsp.query_end == record.query_length)
                                   e2=bool(hsp.query_end >= record.query_length-4)
                                   if (type =='HTH'):
                                        if (s1 and e1):
                                             results[num-1][index]=1
                                   elif (type =='HTT'):
                                        if(s1 and e2):
                                             results[num-1][index]=1
                                   elif (type =='TTH'):
                                        if(s2  and e1):
                                             results[num-1][index]=-1
     print "doneconservation"
     return results

##########################################################
#this is to see if the sequence is conserved in other species, the seq_name_list has some special requiremennt
from Bio.Blast import NCBIXML
from numpy import *
from string import *
def conservation_check3(blast_result_file,conservation_list,seq_name_list):
     results = [[0] * len(conservation_list) for row in range(len(seq_name_list))]
     for record in NCBIXML.parse(open(blast_result_file)) :
          if record.alignments :
               name=record.query
               num=split(record.query,"_")[1]
               type=split(record.query,"_")[2]
               num=int(num,10)
               print num
               for align in record.alignments :
                    for hsp in align.hsps :
                         for index in range(len(conservation_list)):
                              if align.hit_def.find(conservation_list[index]) >-1:
                                   #print "enter"
                                   #b1=bool(hsp.query_start ==1)
                                   #b2=bool(hsp.query_end == record.query_length)
                                   #print "identitiy"
                                   #print hsp.identities
                                   b3=bool(hsp.positives == record.query_length)
                                   b4=bool(hsp.gaps == 0)
                                   if (b3 and b4):
                                        results[num-1][index]=1
     print "doneconservation"
     return results

##########################################################
#this is to see if the sequence is conserved in other species, the seq_name_list has some special requiremennt
from Bio.Blast import NCBIXML
from numpy import *
from string import *
def exclude_check(blast_result_file,exclude_list,seq_name_list):
     results = [[1] * len(exclude_list) for row in range(len(seq_name_list))]
     for record in NCBIXML.parse(open(blast_result_file)) :
          if record.alignments :
          #print >>outfile,record.query,
               name=record.query
               num=split(record.query,"_")[1]
               num=int(num,10)
               #print "num",num
               leng=split(record.query,"_")[2]
               #result=[[0 for x in range(len(exclude_list))] for y in range(128)] #this is a 0/1 matrix to show if the TB type is found in the blast result, there are 4+1=5 genomes to check the exclude and 134 sequences blast results
               #print "result",len(conservation_list)
               for align in record.alignments :
                    for hsp in align.hsps :
                         for index in range(len(exclude_list)):
                              if align.hit_def.find(exclude_list[index]) >-1:
                                   print "into",num
                                   results[num-1][index]=0
     return results
def cal_tm_bond(tri_seq,temp_seq,C_Na,C_Mg,C_Strand):
    tm1 = TmDeltaG.calTm(tri_seq, temp_seq, C_Na, C_Mg,C_Strand, 0.00008)
    tm2 = TmDeltaG.calTm(temp_seq, temp_seq, C_Na, C_Mg,C_Strand, 0.00008)
    bond = str(SecStructures_jf4.SecStructures(SeqRecord(Seq(tri_seq)),SeqRecord(Seq(temp_seq)))).split()[0]
    r=[tm1,tm2,bond]
    return r
###################################################################
def templet_generation_and_prediction(finger_sequence, finger_name,len_minimum, len_maximum, operation_tm, C_Na,C_Mg,C_Strand):
    rank=finger_name.split("_")[1]
    type=finger_name.split("_")[2]
    location=finger_name.split("_")[3]
    length=len(finger_sequence)
    trigger=[]
    template=[]
    maximum=int(len_maximum)
    minimum=int(len_minimum)
    nc=['A','T','G','C']
    xxxx=['']
    if (type=='HTH'):
        tri_gen=finger_sequence[length-9-8:length-9]
        tri_gen_len=len(tri_gen)
        for i in range(minimum,maximum+1):
            if (length-9-i>=9):
                up_trigger=finger_sequence[length-9-i:length-9]
                bottom_trigger=revcomp(finger_sequence[9:i+9])
                for j in range(4):
                    for h in range(4):
                        for k in range(4):
                            for m in range(4):
                                xxxx=nc[j]+nc[h]+nc[k]+nc[m]
                                up_emplate=up_trigger+'GAGTC'+xxxx+up_trigger
                                bottom_template=bottom_trigger+'GAGTC'+xxxx+bottom_trigger
                #print "<tr><td>",up_trigger,"</td><td>",revcomp(up_template),"<td><tr>"
                up_set=cal_tm_bond(up_trigger,revcomp(up_template),1000,0,50)
                up_pred=SeqDep.method_2_prediction(str(up_template).upper())
                bottom_set=cal_tm_bond(up_trigger,revcomp(up_template),C_Na*1000,C_Mg*1000,C_Strand*1000000000)
                bottom_pred=SeqDep.method_2_prediction(str(up_template).upper())
                print up_trigger, revcomp(up_template), up_set[0], configure_finger.min_tri_temp_tm, float(up_set[0]), configure_finger.max_tri_temp_tm, float(up_set[1]), configure_finger.max_temp_tm , int(up_set[2]), configure_finger.max_temp_bonds
                if (float(up_set[0])>configure_finger.min_tri_temp_tm and float(up_set[0])<configure_finger.max_tri_temp_tm and float(up_set[1])<configure_finger.max_temp_tm and int(up_set[2])<configure_finger.max_temp_bonds):
                    todo.Tritemp(finger_id=rank,type=type, start= location, trigger =up_trigger , trig_gen=tri_gen,temp=up_template,tri_length=tri_gen_len ,temp_bayes_class=up_pred[1],temp_pwm_class=up_pred[0][0] ,temp_p90_score=up_pred[0][1],temp_diff_score= up_pred[0][2],tri_temp_tm=up_set[0] ,temp_tm=up_set[1] ,bonds=up_set[2] )
                    elixir.session.commit
                if (float(bottom_set[0])>configure_finger.min_tri_temp_tm and float(bottom_set[0])<configure_finger.max_tri_temp_tm and float(bottom_set[1])<configure_finger.max_temp_tm and int(bottom_set[2])<configure_finger.max_temp_bonds):
                    todo.Tritemp(finger_id=rank,type=type, start= location, trigger =bottom_trigger , trig_gen=tri_gen,temp=bottom_template, tri_length=tri_gen_len ,temp_bayes_class=bottom_pred[1],temp_pwm_class=bottom_pred[0][0] ,temp_p90_score=bottom_pred[0][1],temp_diff_score= bottom_pred[0][2],tri_temp_tm=bottom_set[0] ,temp_tm=bottom_set[1] ,bonds=bottom_set[2] )
                    elixir.session.commit
                
    if (type=='HTT'):
        tri_gen=finger_sequence[length-9-8:length-4]
        len_tri_gen=len(tri_gen)
        #print "HTT enter"
        for i in range(minimum,maximum+1):
            if (length-9-i>=9):
                up_trigger=finger_sequence[length-9-i:length-9]
                xxxx=finger_sequence[length-4:length]
                up_template=up_trigger+'GAGTC'+xxxx+up_trigger
                up_set=cal_tm_bond(up_trigger,revcomp(up_template),C_Na*1000,C_Mg*1000,C_Strand*1000000000)
                up_pred=SeqDep.method_2_prediction(str(up_template).upper())
                if (float(up_set[0])>configure_finger.min_tri_temp_tm and float(up_set[0])<configure_finger.max_tri_temp_tm and float(up_set[1])<configure_finger.max_temp_tm and int(up_set[2])<configure_finger.max_temp_bonds):
                    todo.Tritemp(finger_id=rank,type=type, start= location, trigger =up_trigger , trig_gen=tri_gen,temp=bottom_template,tri_length=tri_gen_len ,temp_bayes_class=up_pred[0],temp_pwm_class=up_pred[1] ,temp_p90_score=up_pred[0][1],temp_diff_score= up_pred[0][2],tri_temp_tm=up_set[0] ,temp_tm=up_set[1] ,bonds=up_set[2] )
                    elixir.session.commit


    if (type=='TTH'):
        for i in range(minimum,maximum+1):
            if (length-9-i>=9):
                bottom_trigger=revcomp(finger_sequence[9:i+9])
                xxxx=revcomp(finger_sequence[0:4])
                bottom_template=bottom_trigger+'GAGTC'+xxxx+bottom_trigger
                bottom_set=cal_tm_bond(bottom_trigger,revcomp(bottom_template),C_Na*1000,C_Mg*1000,C_Strand*1000000000) 
                bottom_pred=SeqDep.method_2_prediction(str(bottom_template).upper())
                if (float(bottom_set[0])>configure_finger.min_tri_temp_tm and float(bottom_set[0])<configure_finger.max_tri_temp_tm and float(bottom_set[1])<configure_finger.max_temp_tm and int(bottom_set[2])<configure_finger.max_temp_bonds):
                    todo.Tritemp(finger_id=rank,type=type, start= location, trigger =bottom_trigger , trig_gen=tri_gen,temp=bottom_template, tri_length=tri_gen_len ,temp_bayes_class=bottom_pred[1],temp_pwm_class=bottom_pred[0][0] ,temp_p90_score=bottom_pred[0][1],temp_diff_score= bottom_pred[0][2],tri_temp_tm=bottom_set[0] ,temp_tm=bottom_set[1] ,bonds=bottom_set[2] )
                    elixir.session.commit


##############################################
# this is for calculate the number of helix and bonds for each helix from .rnaml file which is generated by program UNAfold.pl and ct2rnaml

##################################################
#this is the calculate the multiple result of a list
def mul_m(results,i):
    m=1
    t1=len(results[i])-1
    for j in range(t1):
        m=m*results[i][j]
    return m

####################################################
def get_data(model):
    f=open('taxid9_short.txt').read().split('\n')
    model.setStringList(f)
    #model.setStringList(['above','below','citrus'])
####################################################
def main():
    test=1
    # Initialize database
    # tritemp.initDB()
    # Create two tags
    # Create a few tags and tag them
    #tritemp.Tritemp(pair_id =26 ,trigger =u'ATCGGACTAGCA' , trig_gen=u'ATCGGACTAGCATTAG',temp=u'ATCGGACTAGCAGACTCATCGGACTAGCA' ,tri_length=12 ,temp_bayes_class=u'good' ,temp_pwm_class=u'bad' ,temp_p90_score=8.987 ,temp_diff_score= 9.765,tri_temp_tm=34.2 ,temp_tm=40.1 ,bonds=10 )
    #tritemp.saveData()
    #print "Tasks with l:"
    #print [(t.pair_id,t.trigger,t.trig_gen, t.tri_length, t.temp_p90_score) for t in tritemp.Tritemp.query.filter().all()]

if __name__ == "__main__":
    main()



