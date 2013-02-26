from Bio import SeqIO
from Bio import Entrez
from PyQt4 import QtGui
import elixir
import todo
from xlwt import Workbook
from Bio.Blast import NCBIWWW
import re
import string
from Bio.Blast import NCBIXML
def get_seqid(seq_id):
    try:
        Entrez.email = "jqian@kgi.edu"
        handle = Entrez.efetch(db="nucleotide", rettype="fasta", id=seq_id)
        seq_record = SeqIO.read(handle, "fasta")
        handle.close()
        f=open(seq_id+'.txt','w')
        print >>f, '>'+seq_record.id+'\n'+seq_record.seq
        f.close()
        return seq_id+'.txt'
    except:
        QtGui.QMessageBox.warning(None, "warning", "not an accession number", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
        return ''

def input_content(input_seq,input_file):
    if input_seq!='' and input_file!='':
        QtGui.QMessageBox.warning(None, "warning", "can not both exists", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
        return ''
    else:
        if input_seq=='' and input_file=='':
            QtGui.QMessageBox.warning(None, "warning", "can not both blank", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
            return ''
        else:
            if input_seq!='':
                if re.findall('^>', input_seq):
                    f=open('~seq.fasta','w')
                    print >>f, input_seq
                    f.close()
                    return 'test.fasta'
                else:
                    return get_seqid(input_seq) #input is the accession number of the sequence
            else:
                try:
                    input_filename = open(input_file,'r')
                    return input_filename
                except:
                    QtGui.QMessageBox.warning(None, "warning", "cann't access to the file", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No | QtGui.QMessageBox.Yes)
                    return ''

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

#this is to calculate the gc content of a dna sequence

def GC_content(dna):
    gc  =  (string.count(dna,  'C')  +  string.count(dna,  'c')  +  string.count(dna,  'G')  +  string.count(dna,  'g'))  /  float(len(dna))  *  100
    return gc

###########################################
#this is to do the complementary sequence for dna
def complement(dna):
    t=string.maketrans("AGCTagct",	"TCGAtcga") # to make the translation table
    return translate(dna, t)

################################
#Restriction site occurrences as a list (),find  all  start  positions  of  a  restriction  site
def  restrict(dna,  enz):
    res  =  []
    site = dna.find(enz)
    while (site != -1):
        res.append(site)
        site  =  dna.find(enz,  site  +  1)
    return  res

#####################################
#Get the codon list from a DNA sequence ()
def codons(s,frame=0):
#frame can be 0, 1 or 2
    codons=[]
    end=len(s[frame:]) - (len(s[frame:]) % 3) - 1
    for i in range(frame,end,3):
        codons.append(s[i:i+3])
    return codons

##########################################
#Reverse Complement of DNA,reverse  complement  of  a  DNA  sequence
def  revcomp(dna):
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
def find_finger_seq(seq,list1, list2,filename, maximum,minimum,HTH,HTT,TTH):
#to find the location for all the recognition site and write the trigger to filename
    list1.sort()
    list2.sort()
    rank=1
    outfile=open(filename,'w')
    if (HTT):
        for i in range(len(list1)-1):
            if ((list1[i+1]-list1[i])<=(maximum+9)) and ((list1[i+1]-list1[i])>=(minimum+9)):# this is HTT
                location=list1[i]+1
                length=list1[i+1]-list1[i]+9
                finger=seq[list1[i]:list1[i+1]+9]
                finger_inside=seq[list1[i]+9:list1[i+1]]
                if finger_inside.find('GAGTC')==-1 and finger_inside.find('GACTC')==-1:
                    print >>outfile,">seq_%s_HTT_%s_%s" %(rank,location,length)
                    print >>outfile,finger
                    rank=rank+1
    
    
    if (TTH):
        for i in range(len(list2)-1):
            if ((list2[i+1]-list2[i])<=(maximum+9)) and ((list2[i+1]-list2[i])>=(minimum+9)): # this is TTH
                    location=list2[i]+1-4
                    length=list2[i+1]-list2[i]+9
                    finger=seq[list2[i]-4:list2[i+1]+5]
                    finger_inside=seq[list2[i]+5:list2[i+1]-4]
                    if finger_inside.find('GAGTC')==-1 and finger_inside.find('GACTC')==-1:
                            print >>outfile,">seq_%s_TTH_%s_%s" %(rank,location,length)
                            print >>outfile,finger
                            rank=rank+1
    
    
    if (HTH):
        p=0
        for i in range(len(list1)):
            for j in range(p,len(list2)):
                if list2[j]<list1[i]:
                    p=j
                    continue    
                elif list2[j]-list1[i]>(maximum+13) :
                    break    
                elif ((list2[j]-list1[i])<=(maximum+13)) and ((list2[j]-list1[i])>=(minimum+13)) and (i==len(list1)-1 or list2[j]<list1[i+1]) and (j==0 or list2[j-1]<list1[i]): #this is HTH
                            location=list1[i]+1
                            length=list2[j]-list1[i]+5
                            finger_inside=seq[list1[i]+5:list2[j]]
                            if finger_inside.find('GAGTC')==-1 and finger_inside.find('GACTC')==-1:
                                print >>outfile,">seq_%s_HTH_%s_%s" %(rank,location,length)
                                print >>outfile,seq[list1[i]:list2[j]+5]
                                rank=rank+1
    
    
    outfile.close()

########################
def find_finger(seq,list1, list2):
#to find the location for all the recognition site and write the trigger to file
        pair=[]
        list1.sort()
        list2.sort()
        rank=1
        print  "<table border=3>\r\n"
        print  "<tr>\r\n"
        print  "<td>" ,"rank" ,"</td><td>","location","</td><td>","sequence","</td>"
        print  "</tr>\r\n"
        for i in range(len(list1)-1):
                if ((list1[i+1]-list1[i])<=39 and (list1[i+1]-list1[i])>=17): # this is HTT
                        pair.append(list1[i])
                        pair.append(list1[i+1])
                        #print ">seq_",rank,"_HTT\n",seq[list1[i]:list1[i+1]+9]
                        location=list1[i]+1
                        print  "<tr>"
                        print  "<td>", rank ,"</td>"
                        print  "<td>", location ,"</td><td>"
                        print  "<font color=\"0000FF\">",seq[list1[i]:list1[i]+5]
                        print  "<font color=\"7777777\">",seq[list1[i]+5:list1[i]+9]
                        print  "<font color=\"000000\">",seq[list1[i]+10:list1[i+1]]
                        print  "<font color=\"0000FF\">",seq[list1[i+1]:list1[i+1]+5]
                        print  "<font color=\"7777777\">",seq[list1[i+1]+5:list1[i]+9]
                        print  "</td></tr>"
                        rank=rank+1


        for i in range(len(list2)-1):
                if ((list2[i+1]-list2[i])<=39 and (list2[i+1]-list2[i])>=17): # this is TTH
                        pair.append(list2[i])
                        pair.append(list2[i+1])
                        location=list2[i]+1
                        #print ">seq_",rank,"_TTH\n",seq[list1[i]-4:list1[i+1]+5]
                        print  "<tr>"
                        print  "<td>", rank ,"</td>"
                        print  "<td>", location ,"</td><td>"
                        print  "<font color=\"7777777\">",seq[list2[i]-4:list2[i]]
                        print  "<font color=\"0000FF\">",seq[list2[i]:list2[i]+5]
                        print  "<font color=\"000000\">",seq[list2[i]+5:list2[i+1]-4]
                        print  "<font color=\"7777777\">",seq[list2[i+1]-4:list2[i+1]]
                        print  "<font color=\"0000FF\">",seq[list2[i+1]:list2[i+1]+5]
                        print  "</td></tr>"
                        rank=rank+1


        for i in range(len(list1)):
                for j in range(len(list2)):
                        if (list2[j]-list1[i]<=43) and (list2[j]-list1[i]>=21) and (list2[j]<list1[i+1] or i==len(list2)) and (list2[j-1]<list1[i] or j==0): #this is HTH
                                pair.append(list1[i])
                                pair.append(list2[j])
                                location=list1[i]+1
                                #print ">seq_",rank,"_HTH\n",seq[list1[i]:list1[i+1]+5]
                                print  "<tr>"
                                print  "<td>", rank ,"</td>"
                                print  "<td>", location ,"</td><td>"
                                print  "<font color=\"0000FF\">",seq[list1[i]:list1[i]+5]
                                print  "<font color=\"7777777\">",seq[list1[i]+5:list1[i]+9]
                                print  "<font color=\"000000\">",seq[list1[i]+9:list2[j]]
                                print  "<font color=\"7777777\">",seq[list2[j]-4:list2[j]]
                                print  "<font color=\"0000FF\">",seq[list2[j]:list2[j]+5]
                                print  "</td></tr>"
                                rank=rank+1



        print  "</table>"
        #outfile.close()
        #return pair


########################
#get a sequence by accname from NCBI and write in in fasta format then store it in a file
def get_seq(accname,filename):
    Entrez.email = 'jqian@kgi.edu'
    handle = Entrez.efetch(db='nucleotide',id=accname, rettype='fasta')
    record = SeqIO.read(handle, 'fasta')
    output_file = open(filename,'wb')
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
    typ=[]
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        typ.append(seq_record.id)
    return typ


##############################################
#read in a dictionary between sequence name and rank
def read_diction_name(fasta_file):
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
     #print "doneconservation"
     return results

##########################################################
#this is to see if the sequence is conserved in other species, the seq_name_list has some special requiremennt
def exclude_check_onlycheckname(blast_result_file,exclude_list,seq_name_list):
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

#################################################################
def print_html_find_finger_detail(seq,list1,list2,conservation_result,exclude_result,maximum,minimum,HTH,HTT,TTH):
    list1.sort()
    list2.sort()
    rank=1
    print  "<table border=3>\r\n"
    print  "<tr>\r\n"
    print  "<td>" ,"rank" ,"</td><td>","location","</td><td>","type","</td><td>sequence","</td><td>","conservation_check","</td><td>","exclude_check","</td>"
    print  "</tr>\r\n"
    for i in range(len(list1)-1):
        if ((list1[i+1]-list1[i])<=(int(maximum,10)+9) and (list1[i+1]-list1[i])>=(int(minimum,10)+9)):# this is HTT
            #print ">seq_",rank,"_HTT\n",seq[list1[i]:list1[i+1]+9]
            location=list1[i]+1
            print  "<tr>"
            print  "<td>", rank ,"</td>"
            print  "<td>", location ,"</td><td>"
            print  "HTT</td><td>"
            print  "<font color=\"0000FF\">",seq[list1[i]:list1[i]+5]
            print  "<font color=\"7777777\">",seq[list1[i]+5:list1[i]+9]
            print  "<font color=\"000000\">",seq[list1[i]+10:list1[i+1]]
            print  "<font color=\"0000FF\">",seq[list1[i+1]:list1[i+1]+5]
            print  "<font color=\"7777777\">",seq[list1[i+1]+5:list1[i+1]+9]
            print  "</td>"
            m=1
            for e in range(len(conservation_result[rank][:])):
                m=m*conservation_result[rank][e]

            if m==0:
                print "<td>no</td>"
            else:
                print "<td>ok</td>"

            m=1
            for e in range(len(exclude_result[rank][:])):
                m=m*exclude_result[rank][e]
                #if exclude_result[rank][e]==0:
                    #a=exclude_list[e-1]
            if m==0:
                print "<td>no</td>"
            else:
                print "<td>ok</td>"

            rank=rank+1
            print "</tr>"

    for i in range(len(list2)-1):
        if ((list2[i+1]-list2[i])<=int(maximum,10)+9 and (list2[i+1]-list2[i])>=int(minimum,10)+9): # this is TTH
            location=list2[i]+1
        #print ">seq_",rank,"_TTH\n",seq[list1[i]-4:list1[i+1]+5]
            print  "<tr>"
            print  "<td>", rank ,"</td>"
            print  "<td>", location ,"</td><td>"
            print  "TTH</td><td>"
            print  "<font color=\"7777777\">",seq[list2[i]-4:list2[i]]
            print  "<font color=\"0000FF\">",seq[list2[i]:list2[i]+5]
            print  "<font color=\"000000\">",seq[list2[i]+5:list2[i+1]-4]
            print  "<font color=\"7777777\">",seq[list2[i+1]-4:list2[i+1]]
            print  "<font color=\"0000FF\">",seq[list2[i+1]:list2[i+1]+5]
            print  "</td>"
            m=1
            for e in range(len(conservation_result[rank][:])):
                m=m*conservation_result[rank][e]

            if m==0:
                print "<td>no</td>"
            else:
                print "<td>ok</td>"

            m=1
            for e in range(len(exclude_result[rank][:])):
                m=m*exclude_result[rank][e]
                #if exclude_result[rank][e]==0:
                    #a=exclude_list[e-1]
            if m==0:
                print "<td>no</td>"
            else:
                print "<td>ok</td>"
            print "</tr>"
            rank=rank+1

    for i in range(len(list1)):
        for j in range(len(list2)):
            if (list2[j]-list1[i]<=int(maximum,10)+13) and (list2[j]-list1[i]>=int(minimum,10)+13) and (list2[j]<list1[i+1] or i==len(list2)) and (list2[j-1]<list1[i] or j==0): #this is HTH
                location=list1[i]+1
                #print ">seq_",rank,"_HTH\n",seq[list1[i]:list1[i+1]+5]
                print  "<tr>"
                print  "<td>", rank ,"</td>"
                print  "<td>", location ,"</td><td>"
                print  "HTH</td><td>"
                print  "<font color=\"0000FF\">",seq[list1[i]:list1[i]+5]
                print  "<font color=\"7777777\">",seq[list1[i]+5:list1[i]+9]
                print  "<font color=\"000000\">",seq[list1[i]+9:list2[j]]
                print  "<font color=\"7777777\">",seq[list2[j]-4:list2[j]]
                print  "<font color=\"0000FF\">",seq[list2[j]:list2[j]+5]
                print  "</td>"
                m=1
                for e in range(len(conservation_result[rank][:])):
                    m=m*conservation_result[rank][e]
                if m==0:
                    print "<td>no</td>"
                else:
                    print "<td>ok</td>"

                m=1
                for e in range(len(exclude_result[rank][:])):
                    m=m*exclude_result[rank][e]
                    #if exclude_result[rank][e]==0:
                        #a=exclude_list[e-1]
                if m==0:
                    print "<td>no</td>"
                else:
                    print "<td>ok</td>"
                rank=rank+1
                print "</tr>"

    print  "</table>"

####################################################################

def print_site(seq,type):
    if (type=='HTH'):
        print  "<font color=\"0000FF\">",seq[0:5]
        print  "<font color=\"7777777\">",seq[5:9]
        print  "<font color=\"000000\">",seq[9:(len(seq)-9)]
        print  "<font color=\"7777777\">",seq[(len(seq)-9):(len(seq)-5)]
        print  "<font color=\"0000FF\">",seq[(len(seq)-5):len(seq)]
    elif (type=='TTH'):
        print  "<font color=\"7777777\">",seq[0:4]
        print  "<font color=\"0000FF\">",seq[4:9]
        print  "<font color=\"000000\">",seq[9:(len(seq)-9)]
        print  "<font color=\"7777777\">",seq[(len(seq)-9):(len(seq)-5)]
        print  "<font color=\"0000FF\">",seq[(len(seq)-5):len(seq)]
    elif (type=='HTT'):
        print  "<font color=\"0000FF\">",seq[0:5]
        print  "<font color=\"7777777\">",seq[5:9]
        print  "<font color=\"000000\">",seq[9:(len(seq)-9)]
        print  "<font color=\"0000FF\">",seq[(len(seq)-9):(len(seq)-4)]
        print  "<font color=\"7777777\">",seq[(len(seq)-4):len(seq)]
    else:
        print "wrong type"


###################################################################

##############################################
# this is for calculate the number of helix and bonds for each helix from .rnaml file which is generated by program UNAfold.pl and ct2rnaml
def cal_helix(file):
    input_file=open(file, 'r')
    site=-1
    leng=[]
    count=0
    all_cal_helix=[]
    line=input_file.readline()
    while (line!=''):
        line=input_file.readline()
        if (line.find('seq-data')!=-1):
            count=0
            leng=[]
            #print "new"
        #if (line.find('<model id')!=-1):
            #a=line
        elif (line.find('End of folding 1 for')!=-1):
            #print "end"
            count=count/2
            #print leng
            all_cal_helix.append([count,leng])
            count=0
            leng=[]
            continue
        elif (line.find('helix')!=-1): #count the number of helix
            count=count+1
        elif (line.find('length')!=-1): #count the length of bonds in this helix
            t1=line.find('>')
            t2=line.find('</')
            leng.append(line[t1+1:t2])


    return all_cal_helix
    input_file.close()

########################################################
# this is for calculate the number of helix and bonds for each helix from .rnaml file which is generated by program UNAfold.pl and ct2rnaml,also with the recorded name
def cal_helix_name(file,all_cal_helix):
    input_file=open(file, 'r')
    site=-1
    leng=[]
    count=0
    line=input_file.readline()
    while (line!=''):
        line=input_file.readline()
        if (line.find('seq-data')!=-1):
            count=0
            leng=[]
            #print "new"
        if (line.find('<molecule id')!=-1):
            a=line.find('" type')
            name=line[16:a]
        elif (line.find('End of folding 1 for')!=-1):
            #print "end"
            count=count/2
            #print leng
            if (not(name in all_cal_helix)):
                all_cal_helix[name]=[count,leng]
            count=0
            leng=[]
            name=""
            continue
        elif (line.find('helix')!=-1): #count the number of helix
            count=count+1
        elif (line.find('length')!=-1): #count the length of bonds in this helix
            t1=line.find('>')
            t2=line.find('</')
            leng.append(line[t1+1:t2])


    return all_cal_helix
    input_file.close()

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
    f=open("taxid9_short.txt").read().split('\n')
    model.setStringList(f)
    #model.setStringList(['above','below','citrus'])
####################################################
def submit_to_query(line): #from the input blank to generate the entrez query for NCBIWWW.qblast
    strlist=str(line).split(' OR ')
    txid_line=''
    for valist in strlist:
        txid_num=valist[valist.find('(taxid:')+7:valist.find(')')]
        txid='txid'+txid_num+' [ORGN]'
        if txid_line != '':
            txid_line=txid_line+' OR '+txid
        else:
            txid_line=txid
    return txid_line

####################################################
def include_check(blast_result_filename,include_line,seq_name_list,num):
        strlist=str(include_line).split(' OR ')
        results={}
        for valist in strlist:
                txid_num=valist[valist.find('(taxid:')+7:valist.find(')')]
                blast_result_file= open(blast_result_filename+txid_num,"r")
                found={}
                for record in NCBIXML.parse(blast_result_file):
                        name=record.query
                        if record.alignments :
                                query_len=record.query_letters # the length of your submitted sequence
                                for align in record.alignments :
                                        for hsp in align.hsps :
                                                #print "blast: ",hsp.identities,name,query_len,int(num)
                                                if hsp.identities >= query_len-int(num): # equal or less mistmatch
                                                        if not found.has_key(record.query):
                                                                found[record.query]=1
                                                                if results.has_key(name):
                                                                        temp=results[name]
                                                                        temp.append(txid_num)
                                                                        results[name]=temp
                                                                else:
                                                                        temp=[txid_num]
                                                                        results[name]=temp
                                                        #print name,results[name]
        return (results)
####################################################
def exclude_check(blast_result_file,exclude_line,seq_name_list,perct, GACTC_YES):
        exclude_line2=re.sub('\(taxid:[0-9]+\)','',exclude_line)
        exclude_list=exclude_line2.split(' OR ')
        results = {}
        blast_result_file_handle = open(blast_result_file)
        for record in NCBIXML.parse(blast_result_file_handle) :
                name=record.query # name is the submitted sequence name
                        #                               GACTC_start     GACTC_end
                        #GAGTC.....GAGTC....   HTT      1       end-8
                        #....GACTC......GACTC  TTH      5       end-4
                        #GAGTC..........GACTC  HTH      1       end-4
                if name.split('_')[2] == 'TTH':
                        GACTC_start=5
                else:
                        GACTC_start=1
                if name.split('_')[2] == 'HTT':
                        GACTC_end=int(name.split('_')[4])-8
                else:
                        GACTC_end=int(name.split('_')[4])-4
                found={} #the organisme found to be cross react with in the database
                results[name]='' #begin with not found anything yet
                if record.alignments :
                        query_len=record.query_letters
                        for align in record.alignments :
                                for hsp in align.hsps :
                                        #print "hsp:",hsp.identities, query_len, perct
                                        if hsp.identities >query_len*perct:
                                                if GACTC_YES:
                                                        print "GACTC_YES:",hsp.query_start, GACTC_start, hsp.query_end
                                                        if (hsp.query_start<=GACTC_start and GACTC_start+5<=hsp.query_end) or (hsp.query_start<=GACTC_end and GACTC_end+5<=hsp.query_end):
                                                                if not found.has_key(align.hit_def):
                                                                        print name, results[name]
                                                                        results[name]=results[name]+align.hit_def  # hit_def is the organism found cross react in the NCBI database
                                                                else:
                                                                        results[name]=results[name]+align.hit_def
                                                                        print name, results[name]
        return (results)
        #return(found)


####################################################
def exclude_check716(blast_result_file,exclude_line,seq_name_list,perct, GACTC_YES):
        blast_result_file_handle = open(blast_result_file)
        for record in NCBIXML.parse(blast_result_file_handle) :
                name=record.query # name is the submitted sequence name
                if name.split('_')[2] == 'TTH':
                        GACTC_start=5
                else:
                        GACTC_start=1
                found={} #the organisme found to be cross react with in the database
                results[name]='not found' #begin with not found anything yet
                if record.alignments :
                        query_len=record.query_letters
                        for align in record.alignments :
                                for hsp in align.hsps :
                                        if GACTC_YES:
                                                GACTC_start=int(name[name.find('GACTC')+5:])
                                                if (hsp.query_start<GACTC_start,GACTC_start+5<hsp.query_end):
                                                        if not found.has_key(align.hit_def):
                                                                found[align.hit_def]=1 # hit_def is the organism found cross react in the NCBI database
                                                        if results.has_key(name):
                                                                results[name]='found'
                                                        else:
                                                                results[name]='found'
                                                else:
                                                        if results.has_key(name):
                                                                #results[name]=results[name]+','+align.hit_def
                                                                results[name]='found'
                                                        else:
                                                                #results[name]=align.hit_def
                                                                results[name]='found'
        return (results)
####################################################
def check_blast(input_file,taxid_line,num, GACTC_YES,ex_in,blast_filename,perct):
        query_line=submit_to_query(taxid_line)
        input_file = open(input_file,'r')
        typ=''
        seq_name_list=[]
        #blast_result_file= open(blast_filename,"w")
        for seq_record in SeqIO.parse(input_file, "fasta"):
                if (len(typ)>5000):
                        #result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=query_line,expect=10)
                        #t=result_handle.read()
                        #blast_result_file.write(result_handle.read())
                        #blast_result_file.close()
                        typ=''
                        print "5k done!"
                typ=typ+seq_record.format('fasta')
                seq_name_list.append(seq_record.id)
        #blast_result_file.close()
        #blast_result_file= open(blast_filename,"a")
        #result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=query_line,expect=10)
        #t2=result_handle.read()
        #blast_result_file.write(t2)
        #blast_result_file.close()
        input_file.close()
        print "blast job done!"
        if ex_in=='ex':
            (check_blast_results)=exclude_check(blast_filename,taxid_line,seq_name_list,perct,GACTC_YES)
        elif ex_in=='in':
            (check_blast_results)=include_check(blast_filename,taxid_line,seq_name_list,num)
            print check_blast_results
        return (check_blast_results)
#################################################################33
def check_blast_716(input_file,taxid_line,num, GACTC_YES,ex_in,blast_filename,evalue):
        query_line=submit_to_query(taxid_line)
        input_file = open(input_file,'r')
        typ=''
        seq_name_list=[]
        blast_result_file= open(blast_filename,"w")
        for seq_record in SeqIO.parse(input_file, "fasta"):
                if (len(typ)>5000):
                        result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=query_line,expect=evalue)
                        blast_result_file.write(result_handle.read())
                        typ=''
                        print "5k done!"
                typ=typ+seq_record.format('fasta')
                seq_name_list.append(seq_record.id)
        blast_result_file.close()
        blast_result_file= open(blast_filename,"a")
        result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=query_line,expect=evalue)
        blast_result_file.write(result_handle.read())
        blast_result_file.close()
        input_file.close()
        print "blast job done!"
        if ex_in=='ex':
            (check_blast_results)=exclude_check(blast_filename,taxid_line,seq_name_list,GACTC_YES)
        elif ex_in=='in':
            (check_blast_results)=include_check(blast_filename,taxid_line,seq_name_list)
        return (check_blast_results)
####################################################
def check_blast_ex(input_filename,taxid_line,num, GACTC_YES,blast_filename,perct):
        query_line=submit_to_query(taxid_line)
        #print "ex", query_line
        input_file = open(input_filename,'r')
        typ=''
        blast_result_file= open(blast_filename,"w")
        for seq_record in SeqIO.parse(input_file, "fasta"):
                if (len(typ)>300):
                        result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=query_line,expect=10)
                        t=result_handle.read()
                        blast_result_file.write(t)
                        typ=''
                        #print "300 done!"
                typ=typ+seq_record.format('fasta')
        blast_result_file.close()
        blast_result_file= open(blast_filename,"a")
        result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=query_line,expect=10)
        t2=result_handle.read()
        blast_result_file.write(t2)
        blast_result_file.close()
        input_file.close()
        print "blast job done!"
##      if ex_in=='ex':
##            (check_blast_results)=exclude_check(blast_filename,taxid_line,seq_name_list,perct,GACTC_YES)
##    elif ex_in=='in':
##            (check_blast_results)=include_check(blast_filename,taxid_line,seq_name_list,num)
##      return (check_blast_results)
####################################################
def check_blast_in(input_filename,taxid_line,num, GACTC_YES,blast_filename,perct):
        #print "in", taxid_line
        strlist=str(taxid_line).split(' OR ')
        for valist in strlist:
                txid_num=valist[valist.find('(taxid:')+7:valist.find(')')]
                blast_result_file= open(blast_filename+txid_num,"w")
                #print "in",blast_result_file
                txid='txid'+txid_num+' [ORGN]'
                #blast_result_file.write('<TAXID>'+txid_num+'</TAXID>'+'\n')
                typ=''
                input_file = open(input_filename,"r")
                for seq_record in SeqIO.parse(input_file, "fasta"):
                        if (len(typ)>200):
                                result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=txid,expect=10)
                                t=result_handle.read()
                                blast_result_file.write(t)
                                typ=''
                                #print "200 done!"
                        typ=typ+seq_record.format('fasta')
                        #print "wating", typ,"finished waiting","\n\n"
                if (len(typ)>0):
                        #print "working on the leftover"
                        result_handle = NCBIWWW.qblast("blastn", "nr", typ,word_size=13,hitlist_size=100,entrez_query=txid,expect=10)
                        t2=result_handle.read()
                        #print typ
                        blast_result_file.write(t2)
                input_file.close()
                blast_result_file.close()
        print "blast job done!"
##    (check_blast_results)=include_check(blast_filename,taxid_line,seq_name_list,num)
##      return (check_blast_results)
####################################################
def combine_check_exinclude(input_file,include_line,exclude_line,num,perct,seq_name_list):
        GACTC_YES=1
        blast_filename='blast_result.xml'
        if exclude_line=='':
                check_blast_results_ex={}
        else:
                check_blast_ex(input_file,exclude_line,num, GACTC_YES,blast_filename,perct)
                (check_blast_results_ex)=exclude_check(blast_filename,exclude_line,seq_name_list,perct,GACTC_YES) #if len(results[name])!=0, then good
        if include_line=='':
                check_blast_results_in={}
        else:
                check_blast_in(input_file,include_line,num, GACTC_YES,blast_filename,perct)
                (check_blast_results_in)=include_check(blast_filename,include_line,seq_name_list,num) #if len(results[name])!=0, then good
        #print "ex:,", check_blast_results_ex
        #print "in:",check_blast_results_in
        return (check_blast_results_in,check_blast_results_ex)
####################################################
def main_page(target_file, include_line, exclude_line,HTH,HTT,TTH,maximum, minimum,num,perct):
        QtGui.QMessageBox.warning(None, "waiting", "It may takes a while,until another window pops out!", QtGui.QMessageBox.Yes)
        for cur_record_o in SeqIO.parse(target_file, "fasta"):
                input_seq=cur_record_o.seq
        list1=restrict(input_seq,  'GAGTC')
        list1.sort()
        list2=restrict(input_seq,  'GACTC')
        list2.sort()
        find_finger_seq(input_seq,list1,list2,"~seq.fasta",maximum,minimum,HTH,HTT,TTH)
        seq_name_list=[]
        for seq_record in SeqIO.parse("~seq.fasta", "fasta"):
                seq_name_list.append(seq_record.id)
        (in_blast_result,ex_blast_result)=combine_check_exinclude('~seq.fasta',str(include_line),str(exclude_line),int(num),float(perct)/100,seq_name_list)
        for cur_record_o in SeqIO.parse('~seq.fasta','fasta'):
                list=cur_record_o.id.split('_')
                if not in_blast_result.has_key(cur_record_o.id):
                        in_blast_result[cur_record_o.id]=''
                if not ex_blast_result.has_key(cur_record_o.id):
                        ex_blast_result[cur_record_o.id]=''
                todo.Task(finger_type=list[2], start=list[3], length=list[4], sequence=str(cur_record_o.seq), include=len(in_blast_result[cur_record_o.id]), exclude=len(ex_blast_result[cur_record_o.id]), done=False)
        elixir.session.commit
def main_page_test(target_file, include_line, exclude_line,HTH,HTT,TTH,maximum, minimum,num,perct):
        for cur_record_o in SeqIO.parse(target_file, "fasta"):
                input_seq=cur_record_o.seq
        list1=restrict(input_seq,  'GAGTC')
        list1.sort()
        list2=restrict(input_seq,  'GACTC')
        list2.sort()
        find_finger_seq(input_seq,list1,list2,"~seq.fasta",maximum,minimum,HTH,HTT,TTH)
def write_xls(filename,selfobj, finger, tritemp):
    book = Workbook()
    sheet1 = book.add_sheet('Parameter',cell_overwrite_ok=True)
    sheet2 = book.add_sheet('All Fingerprintg site',cell_overwrite_ok=True)
    sheet3 = book.add_sheet('Generated Tirgger templates',cell_overwrite_ok=True)
    sheet4 = book.add_sheet('Explanations',cell_overwrite_ok=True)
    row11 = sheet1.row(0)
    row11.write(0,'Input_file')
    row11.write(1,'Include_list')
    row11.write(2,'Exclude_list')
    row11.write(3,'H-H')
    row11.write(4,'H-T')
    row11.write(5,'T-H')
    row11.write(6,'Maximum trigger length')
    row11.write(7,'Minimum trigger length')
    row11.write(8,'Maximum mismatch with include list')
    row11.write(9,'Minimum coverage with exclude list')

    row21 = sheet2.row(0)
    row21.write(0,'Selected')
    row21.write(1,'ID')
    row21.write(2,'Type')
    row21.write(3,'Start')
    row21.write(4,'Length')
    row21.write(5,'Sequence')
    row21.write(6,'Include')
    row21.write(7,'Exclude')


    row31 = sheet3.row(0)
    row31.write(0,'Pair_ID')
    row31.write(1,'Finger_ID')
    row31.write(2,'Type')
    row31.write(3,'Start')
    row31.write(4,'Trigger_generated')
    row31.write(5,'Trigger')
    row31.write(6,'Template')
    row31.write(7,'Bayes_class')
    row31.write(8,'Pwm_class')
    row31.write(9,'P90_score')
    row31.write(10,'Diff_score')
    row31.write(11,'Tri-temp Tm')
    row31.write(12,'Temp_temp_tm')
    row31.write(13,'Bonds')

    row41 = sheet4.row(0)
    row41.write(0,'Pair_ID')
    row41.write(1,'Finger_ID')
    row41.write(2,'Type')
    row41.write(3,'Start')
    row41.write(4,'Length')
    row41.write(5,'Sequence')
    row41.write(6,'Include')
    row41.write(7,'Exclude')
    row41.write(8,'Trigger_generated')
    row41.write(9,'Trigger')
    row41.write(10,'Template')
    row41.write(11,'Bayes_class')
    row41.write(12,'Pwm_class')
    row41.write(13,'P90_score')
    row41.write(14,'Diff_score')
    row41.write(15,'Tri-temp Tm')
    row41.write(16,'Temp_temp_tm')
    row41.write(17,'Bonds')

    row42 = sheet4.row(1)
    row42.write(0,'An Unique Identification number for this trigger and template combination')
    row42.write(1,'An Unique Identification number for the fingerprinting site')
    row42.write(2,'This fingerprinting site is Head to Head (HTH), or Tail to Head(TTH) or Head to Tail(HTT)')
    row42.write(3,'Where does this fingerprinting site locate on user defined sequence')
    row42.write(4,'Length of the finterprinting site')
    row42.write(5,'The sequence of this finerprinting site')
    row42.write(6,'How many species in the include list were found to be matched with the user defined sequence')
    row42.write(7,'How many pairs of matched fragment were found between the user defined sequence and exclude list species')
    row42.write(8,'The trigger generated by the nicking enzyme')
    row42.write(9,'Designed trigger')
    row42.write(10,'Template')
    row42.write(11,'Template performance predicted by bayes method')
    row42.write(12,'Template performance predicted by position weight matrix')
    row42.write(13,'Predicted P90 score for the template, lower P90 score may indicates faster amplification speed')
    row42.write(14,'Predicted Diff score for the template, higher Diff score may indicates longer separation time of the negative and positive sample')
    row42.write(15,'Melting temperature for the designed trigger sequence and template sequence')
    row42.write(16,'Melting temperature for the designed template sequence and template sequence')
    row42.write(17,'Number of predicted secondary bonds between ')



    i1=1
    row1=sheet1.row(i1)
    row1.write(0,str(selfobj.ui.lineEdit.text()))
    row1.write(1,str(selfobj.ui.lineEdit_5.text()))
    row1.write(2,str(selfobj.ui.lineEdit_15.text()))
    row1.write(3,str(selfobj.ui.checkBox.isChecked()))
    row1.write(4,str(selfobj.ui.checkBox_2.isChecked()))
    row1.write(5,str(selfobj.ui.checkBox_3.isChecked()))
    row1.write(6,int(selfobj.ui.lineEdit_10.text()))
    row1.write(7,int(selfobj.ui.lineEdit_9.text()))
    row1.write(8,int(selfobj.ui.lineEdit_2.text()))
    row1.write(9,int(selfobj.ui.lineEdit_3.text()))

    i2=1
    for fr in finger:
                row2=sheet2.row(i2)
                row2.write(0,str(fr.done))
                row2.write(1,str(fr.finger_id))
                row2.write(2,str(fr.finger_type))
                row2.write(3,str(fr.start))
                row2.write(4,str(fr.length))
                row2.write(5,str(fr.sequence))
                row2.write(6,str(fr.include))
                row2.write(7,str(fr.exclude))
                i2=i2+1
    i3=1
    for tt in tritemp:
                row3=sheet3.row(i3)
                row3.write(0,str(tt.pair_id))
                row3.write(1,str(tt.finger_id))
                row3.write(2,str(tt.type))
                row3.write(3,str(tt.start))
                row3.write(4,str(tt.trig_gen))
                row3.write(5,str(tt.trigger))
                row3.write(6,str(tt.temp))
                row3.write(7,str(tt.temp_bayes_class))
                row3.write(8,str(tt.temp_pwm_class))
                row3.write(9,str(tt.temp_p90_score))
                row3.write(10,str(tt.temp_diff_score))
                row3.write(12,str(tt.tri_temp_tm))
                row3.write(12,str(tt.temp_tm))
                row3.write(13,str(tt.bonds))
                i3=i3+1
    book.save(filename)

def write_xls_2(filename,selfobj, finger, tritemp):
    book = Workbook()
    sheet1 = book.add_sheet('Parameter',cell_overwrite_ok=True)
    sheet2 = book.add_sheet('All Fingerprintg site',cell_overwrite_ok=True)
    sheet3 = book.add_sheet('Explanation',cell_overwrite_ok=True)

    row11 = sheet1.row(0)
    row11.write(0,'Input_file')
    row11.write(1,'Include_list')
    row11.write(2,'Exclude_list')
    row11.write(3,'H-H')
    row11.write(4,'H-T')
    row11.write(5,'T-H')
    row11.write(6,'Maximum trigger length')
    row11.write(7,'Minimum trigger length')
    row11.write(8,'Maximum mismatch with include list')
    row11.write(9,'Minimum coverage with exclude list')

    row21 = sheet2.row(0)
    row21.write(0,'Selected')
    row21.write(1,'ID')
    row21.write(2,'Type')
    row21.write(3,'Start')
    row21.write(4,'Length')
    row21.write(5,'Sequence')
    row21.write(6,'Include')
    row21.write(7,'Exclude')

    row31 = sheet3.row(0)
    row31.write(0,'Selected')
    row31.write(1,'ID')
    row31.write(2,'Type')
    row31.write(3,'Start')
    row31.write(4,'Length')
    row31.write(5,'Sequence')
    row31.write(6,'Include')
    row31.write(7,'Exclude')


    row32 = sheet3.row(1)
    row32.write(0,'Whether this Fingerprintig site is selected by the user')
    row32.write(1,'An Unique Identification number for the fingerprinting site')
    row32.write(2,'This fingerprinting site is Head to Head (HTH), or Tail to Head(TTH) or Head to Tail(HTT)')
    row32.write(3,'Where does this fingerprinting site locate on user defined sequence')
    row32.write(4,'Length of the finterprinting site')
    row32.write(5,'The sequence of this finerprinting site')
    row32.write(6,'How many species in the include list were found to be matched with the user defined sequence')
    row32.write(7,'How many pairs of matched fragment were found between the user defined sequence and exclude list species')


    i1=1
    row1=sheet1.row(i1)
    row1.write(0,str(selfobj.ui.lineEdit.text()))
    row1.write(1,str(selfobj.ui.lineEdit_5.text()))
    row1.write(2,str(selfobj.ui.lineEdit_15.text()))
    row1.write(3,str(selfobj.ui.checkBox.isChecked()))
    row1.write(4,str(selfobj.ui.checkBox_2.isChecked()))
    row1.write(5,str(selfobj.ui.checkBox_3.isChecked()))
    row1.write(6,int(selfobj.ui.lineEdit_10.text()))
    row1.write(7,int(selfobj.ui.lineEdit_9.text()))
    row1.write(8,int(selfobj.ui.lineEdit_2.text()))
    row1.write(9,int(selfobj.ui.lineEdit_3.text()))

    i2=1
    for fr in finger:
                row2=sheet2.row(i2)
                row2.write(0,str(fr.done))
                row2.write(1,str(fr.finger_id))
                row2.write(2,str(fr.finger_type))
                row2.write(3,str(fr.start))
                row2.write(4,str(fr.length))
                row2.write(5,str(fr.sequence))
                row2.write(6,str(fr.include))
                row2.write(7,str(fr.exclude))
                i2=i2+1
    book.save(filename)

if __name__ == "__main__":    
    main_page_test('HSV2.fasta','','',1,1,1,30,8,4,66)
