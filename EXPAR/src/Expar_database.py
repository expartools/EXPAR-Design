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
import SecStructures_jf7
from Bio.SeqRecord import SeqRecord
from Bio import Entrez, SeqIO
from Bio.Seq import Seq
import subprocess
import os
import math
#currentPath = os.path.realpath(os.path.normpath(os.path.dirname(__file__)))
#os.environ["PATH"] = os.environ["PATH"]+";"+currentPath+"\\UNAFold\\bin"
#os.environ["PATH"] = os.environ["PATH"]+";.\\UNAFold\\share\\unafold"
#os.sys.path.append(".\\UNAFold\\bin")
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

def cal_tm_bond(tri_seq,temp_seq,C_Na,C_Mg,C_Strand):
    tm1 = TmDeltaG.calTm(tri_seq, temp_seq, C_Na, C_Mg,C_Strand, 0.00008)
    tm2 = TmDeltaG.calTm(temp_seq, temp_seq, C_Na, C_Mg,C_Strand, 0.00008)
    bond = str(SecStructures_jf7.SecStructures(SeqRecord(Seq(tri_seq)),SeqRecord(Seq(temp_seq)))).split()[0]
    r=[tm1,tm2,bond]
    return r
###################################################################
def templet_generation_and_prediction(finger_sequence, finger_name,len_minimum, len_maximum, operation_tm, C_Na,C_Mg,C_Strand):
	if template_trigger_gengerate(finger_sequence,finger_name,len_maximum,len_minimum)==1: #gengerate trigger, template successfully and stored in ~tri.txt, ~tem.txt
		dic_therm=unafold_pred('~tem.txt',"~tri.txt",operation_tm, C_Na,C_Mg,C_Strand)
		dic_jf={}
		dic_tri={}
		dic_temp={}
		dic_trigen={}
		for cur_record_o in SeqIO.parse(open('~tem.txt'), "fasta") :
			dic_temp[str(cur_record_o.id)]=str(cur_record_o.seq)
		for cur_record_o in SeqIO.parse(open('~tri.txt'), "fasta") :
			dic_tri[str(cur_record_o.id)]=str(cur_record_o.seq)
		for cur_record_o in SeqIO.parse(open('~tri_gen.txt'), "fasta") :
			dic_trigen[str(cur_record_o.id)]=str(cur_record_o.seq)
		for cur_record_o in SeqIO.parse(open('~tem.txt'), "fasta") :
			cur_seq=str(cur_record_o.seq).upper() #all capitalized
			dic_jf[cur_record_o.id]=SeqDep.method_2_prediction(cur_seq)
		for key,value in sorted(dic_therm. iteritems(), key=lambda(k,v):(v,k)):
			valuelist=key.split('_')
			part='_'.join(valuelist[0:5])
			part2=part+'_'
			triname=part2+'tri_'
			triname=triname+valuelist[6]
			trigenname=part2+'tri_gen'
			if (float(dic_therm[key][0])>configure_finger.min_tri_temp_tm and float(dic_therm[key][0])<configure_finger.max_tri_temp_tm and float(dic_therm[key][1])<configure_finger.max_temp_tm and int(dic_therm[key][2])<configure_finger.max_temp_bonds):
			#if(float(dic_therm[key][0])>configure_finger.min_tri_temp_tm and float(dic_therm[key][0])<configure_finger.max_tri_temp_tm and float(dic_therm[key][1])<configure_finger.max_temp_tm):
				todo.Tritemp(finger_id=valuelist[1],type=valuelist[2], start= valuelist[3], trigger =dic_tri[triname] , trig_gen=dic_trigen[trigenname],temp=dic_temp[key],tri_length=len(dic_tri[triname]) ,temp_bayes_class=dic_jf[key][1][0],temp_pwm_class=dic_jf[key][0][0] ,temp_p90_score=dic_jf[key][0][1],temp_diff_score= dic_jf[key][0][2],tri_temp_tm=dic_therm[key][0] ,temp_tm=dic_therm[key][1] ,bonds=dic_therm[key][2] )
			elixir.session.commit
##############################################
# this is for calculate the number of helix and melting temperature
def template_trigger_gengerate(sequence,name,maximum=16,minimum=8):
#print "entered"
	tri_file= open("~tri.txt","w")
	tri_gen_file= open("~tri_gen.txt","w")
	tem_file= open("~tem.txt","w")
	rank=name.split("_")[1]
	type=name.split("_")[2]
	location=name.split("_")[3]
	length=len(sequence)
	trigger=[]
	template=[]
	bondadd=0
	num=0
	maximum=int(maximum)
	minimum=int(minimum)
	nc=['A','T','G','C']
	xxxx=['']
	countadd=0
	if (type=='HTH'):
		print >>tri_gen_file, ">"+name+"_tri_gen\n"+sequence[length-9-8:length-9]
		for i in range(minimum,maximum+1):
			if (length-9-i>=9):
				up_trigger=sequence[length-9-i:length-9]
				bottom_trigger=revcomp(sequence[9:i+9])
				for j in range(4):
					for h in range(4):
						for k in range(4):
							for m in range(4):
								xxxx=nc[j]+nc[h]+nc[k]+nc[m]
								up_template=up_trigger+'GAGTC'+xxxx+up_trigger
								bottom_template=bottom_trigger+'GAGTC'+xxxx+bottom_trigger
				#print "<tr><td>",up_trigger,"</td><td>",revcomp(up_template),"<td><tr>"
				trigger.append(up_trigger)
				template.append(revcomp(up_template))
				trigger.append(bottom_trigger)
				template.append(revcomp(bottom_template))
				print >>tri_file, ">"+name+"_tri_"+str(countadd)+"\n"+up_trigger
				print >>tem_file, ">"+name+"_tem_"+str(countadd)+"\n"+revcomp(up_template)
				countadd=countadd+1
				print >>tri_file, ">"+name+"_tri_"+str(countadd)+"\n"+bottom_trigger
				print >>tem_file, ">"+name+"_tem_"+str(countadd)+"\n"+revcomp(bottom_template)
				countadd=countadd+1
	if (type=='HTT'):
		print >>tri_gen_file, ">"+name+"_tri_gen\n"+sequence[length-9-8:length-4]
		#print "HTT enter"
		for i in range(minimum,maximum+1):
			if (length-9-i>=9):
				up_trigger=sequence[length-9-i:length-9]
				xxxx=sequence[length-4:length]
				up_template=up_trigger+'GAGTC'+xxxx+up_trigger
		trigger.append(up_trigger)
		#print "up_template",up_template,'revc',revcomp(up_template)
		template.append(revcomp(up_template))
		print >>tri_file, ">"+name+"_tri_"+str(countadd)+"\n"+up_trigger
		print >>tem_file, ">"+name+"_tem_"+str(countadd)+"\n"+revcomp(up_template)
		countadd=countadd+1


	if (type=='TTH'):
		print >>tri_gen_file, ">"+name+"_tri_gen\n"+sequence[4:17]
		for i in range(minimum,maximum+1):
			if (length-9-i>=9):
				bottom_trigger=revcomp(sequence[9:i+9])
				xxxx=revcomp(sequence[0:4])
				bottom_template=bottom_trigger+'GAGTC'+xxxx+bottom_trigger
				trigger.append(bottom_trigger)
				template.append(revcomp(bottom_template))
				print >>tri_file, ">"+name+"_tri_"+str(countadd)+"\n"+bottom_trigger
				print >>tem_file, ">"+name+"_tem_"+str(countadd)+"\n"+revcomp(bottom_template)
				countadd=countadd+1
	tri_file.close()
	tem_file.close()
	tri_gen_file.close()
	return 1
####################################################################################################
def unafold_pred(temp_filename,trig_filename, operation_tm, C_Na,C_Mg,C_Strand):
    currentPath = os.path.realpath(os.path.normpath(os.path.dirname(__file__)))
    os.chdir(os.path.realpath(os.path.normpath(os.path.dirname(__file__)))+"\\UNAFold\\bin")
    os.popen('del ~*')
    os.popen('del ~*')  
    os.popen('copy '+currentPath+'\\'+trig_filename+' '+currentPath+'.\\UNAFold\\bin')
    os.popen('copy '+currentPath+'\\'+temp_filename+' '+currentPath+'.\\UNAFold\\bin') 
    dic_therm={}
    name=[]
    f=open(temp_filename,'r')        
    num=0
    bond={}
    for cur_record in SeqIO.parse(f, "fasta") :
        name.append(cur_record.id)
        #print str(SecStructures_jf7.SecStructures(cur_record,cur_record)).split()[0]
        bond[cur_record.id]=str(SecStructures_jf7.SecStructures(cur_record,cur_record)).split()[0]    
    f.close()  
    f=open(trig_filename+temp_filename+"DGHS.txt",'w')
    nc="melt.pl -n DNA -t "+str(operation_tm)+" -N "+str(C_Na)+" -M "+str(C_Mg)+" -C "+str(C_Strand)+" "+trig_filename+" "+temp_filename
    p = os.popen(nc)
    print >>f, p.read()
    f.close()
    f=open(temp_filename+temp_filename+"DGHS_2.txt",'w')
    nc2="melt.pl -n DNA -t "+str(operation_tm)+" -N "+str(C_Na)+" -M "+str(C_Mg)+" -C "+str(C_Strand)+" "+temp_filename+" "+temp_filename
    p = os.popen(nc2)
    print >>f, p.read()
    f.close()
    infile = open(trig_filename+temp_filename+"DGHS.txt")
    infile2 = open(temp_filename+temp_filename+"DGHS_2.txt")
    line = infile.readline()
    line2 = infile2.readline()
    i=0
    while 1:
    	line = infile.readline().rstrip()
    	line2 = infile2.readline().rstrip()
    	if line.find('Calculating')==-1 and line2.find('Calculating')==-1 and line.find('dG')==-1 and line2.find('dG')==-1 and line !='' and line2!='':
    	    line=line.replace('-','')
    	    ener=line.split("\t")
    	    ener2=line2.split("\t") #template-template TM
    	    TM=float(ener[1])/((float(ener[2])-(1.987*math.log(C_Strand)))/1000)-273.15 #trigger-template TM
            dic_therm[name[i]]=[TM,float(ener2[3]),bond[name[i]]] # Trigger-template TM, template-template TM,bonds
            num=num+1
            i=i+1
    	if line == "":
    	    break
    infile.close()
    infile2.close()
    os.popen('del ~*')
    os.popen('del ~*')
    print "finished one fingerprint site!"
    os.chdir(currentPath)    
    return dic_therm

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
    #print os.environ["PATH"]
    #currentPath = os.path.realpath(os.path.normpath(os.path.dirname(__file__)))
    #print currentPath
    #subprocess.Popen("melt.pl -n DNA -t 55 -N 1 -M 0 -C 0.0000005 ..\\..\\~tri.txt ..\\..\\~tem.txt", env=my_env)
    #subprocess.Popen("melt.pl -n DNA -t 55 -N 1 -M 0 -C 0.0000005 ..\..\~tri.txt ..\..\~tem.txt", env=my_env)
    currentPath = os.path.realpath(os.path.normpath(os.path.dirname(__file__)))
    #os.system('move ~tri.txt .\\UNAFold\\bin')
    #os.system('move ~tem.txt .\\UNAFold\\bin')
    os.chdir(os.path.realpath(os.path.normpath(os.path.dirname(__file__)))+"\\UNAFold\\bin")
    #print os.path.realpath(os.path.normpath(os.path.dirname(__file__)))
    #os.system('dir')
    nc="melt.pl -n DNA -t 55 -N 1 -M 0 -C 0.0000005 ~tri.txt ~tem.txt"
    p = os.popen(nc)
    f=open('test.txt','w')
    print >>f, p.read()
    f.close()
    
    #os.system('melt.pl -n DNA -t 55 -N 1 -M 0 -C 0.0000005 ..\\..\\src\\~tri.txt ..\\..\\src\\UNAFold\bin\\~tem.txt')
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



