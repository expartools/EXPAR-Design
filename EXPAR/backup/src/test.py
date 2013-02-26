MeltingTemp_jf.calDeltaG(s, t, mono_conc=mono_conc, Mg_conc=Mg_conc, dntp_conc=dntp_conc)


def template_trigger_gengerate(sequence,name,maximum=16,minimum=8):
#this is to generate a template by the finger printting site on the genome
#print "entered"
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
        tri_gen=sequence[length-9-8:length-9]
        tri_gen_len=len(tri_gen)
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
                up_set=cal_tm_bond(up_trigger,revcomp(up_template))
                up_pred=SeqDep.method_2_prediction(cur_seq=str(up_template).upper())
                bottom_set=cal_tm_bond(up_trigger,revcomp(up_template))
                bottom_pred=SeqDep.method_2_prediction(cur_seq=str(up_template).upper())
                if (float(up_set[0])>configure_finger.min_tri_temp_tm and float(up_set[0])<configure_finger.max_tri_temp_tm and float(up_set[1])<configure_finger.max_temp_tm and float(up_set[1])>configure_finger.min_temp_tm and int(up_set[2])<configure_finger.max_temp_bonds):
                    todo.Tritemp(finger_id=rank,type=type, start= location, trigger =up_trigger , trig_gen=tri_gen,temp=up_template,tri_length=tri_gen_len ,temp_bayes_class=up_pred[1],temp_pwm_class=up_pred[0][0] ,temp_p90_score=up_pred[0][1],temp_diff_score= up_pred[0][2],tri_temp_tm=up_set[0] ,temp_tm=up_set[1] ,bonds=up_set[2] )
                    elixir.session.commit
                if (float(bottom_set[0])>configure_finger.min_tri_temp_tm and float(bottom_set[0])<configure_finger.max_tri_temp_tm and float(bottom_set[1])<configure_finger.max_temp_tm and float(bottom_set[1])>configure_finger.min_temp_tm and int(bottom_set[2])<configure_finger.max_temp_bonds):
                    todo.Tritemp(finger_id=rank,type=type, start= location, trigger =bottom_trigger , trig_gen=tri_gen,temp=bottom_template, tri_length=tri_gen_len ,temp_bayes_class=bottom_pred[1],temp_pwm_class=bottom_pred[0][0] ,temp_p90_score=bottom_pred[0][1],temp_diff_score= bottom_pred[0][2],tri_temp_tm=bottom_set[0] ,temp_tm=bottom_set[1] ,bonds=bottom_set[2] )
                    elixir.session.commit
                
    if (type=='HTT'):
        tri_gen=sequence[length-9-8:length-4]
        len_tri_gen=len(tri_gen)
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
        up_set=cal_tm_bond(up_trigger,revcomp(up_template))
        up_pred=SeqDep.method_2_prediction(cur_seq=str(up_template).upper())
        if (float(up_set[0])>configure_finger.min_tri_temp_tm and float(up_set[0])<configure_finger.max_tri_temp_tm and float(up_set[1])<configure_finger.max_temp_tm and float(up_set[1])>configure_finger.min_temp_tm and int(up_set[2])<configure_finger.max_temp_bonds):
            todo.Tritemp(finger_id=rank,type=type, start= location, trigger =up_trigger , trig_gen=tri_gen,temp=bottom_template,tri_length=tri_gen_len ,temp_bayes_class=up_pred[0],temp_pwm_class=up_pred[1] ,temp_p90_score=up_pred[0][1],temp_diff_score= up_pred[0][2],tri_temp_tm=up_set[0] ,temp_tm=up_set[1] ,bonds=up_set[2] )
            elixir.session.commit


    if (type=='TTH'):
        print >>tri_gen_file, ">"+name+"_tri_gen\n"+sequence[4:17]
        for i in range(minimum,maximum+1):
            if (length-9-i>=9):
                bottom_trigger=revcomp(sequence[9:i+9])
                xxxx=revcomp(sequence[0:4])
                bottom_template=bottom_trigger+'GAGTC'+xxxx+bottom_trigger
                trigger.append(bottom_trigger)
                template.append(revcomp(bottom_template))
                bottom_set=cal_tm_bond(up_trigger,revcomp(up_template))
                bottom_pred=SeqDep.method_2_prediction(cur_seq=str(up_template).upper())
                if (float(bottom_set[0])>configure_finger.min_tri_temp_tm and float(bottom_set[0])<configure_finger.max_tri_temp_tm and float(bottom_set[1])<configure_finger.max_temp_tm and float(bottom_set[1])>configure_finger.min_temp_tm and int(bottom_set[2])<configure_finger.max_temp_bonds):
                    todo.Tritemp(finger_id=rank,type=type, start= location, trigger =bottom_trigger , trig_gen=tri_gen,temp=bottom_template, tri_length=tri_gen_len ,temp_bayes_class=bottom_pred[1],temp_pwm_class=bottom_pred[0][0] ,temp_p90_score=bottom_pred[0][1],temp_diff_score= bottom_pred[0][2],tri_temp_tm=bottom_set[0] ,temp_tm=bottom_set[1] ,bonds=bottom_set[2] )
                    elixir.session.commit