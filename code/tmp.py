#set fdm=marker

#################################################################
#Creating Paolas Special tables {{{

#################################################################
#START creating the MRM Atlas tables {{{


#first link the MRM atlas table to the yeast genome
"""
drop table MRMPepLink;
create table MRMPepLink (
 peptide_key int,
 mrm_key int,
 charge int
);
alter table MRMPepLink add index(peptide_key);
"""


c.execute( 'select id,sequence from hroest.MRMAtlas_yeast_qqq_201002_q2')
mrm = c.fetchall()
c.execute( 'select id,sequence from ddb.peptide where experiment_key = 3131')
pep = c.fetchall()



#create a mapping of sequence to peptide_id
pepd = {}
for p in pep:
    pepd[p[1]] = p[0]


#TODO filter out things like PEPT[160]IDE
import re 
res = []
notfound_dict = {}
nontr = 0
other = 0
for m in mrm:
    try:
        spl = m[1].split('/')
        sequence = re.sub('[^A-Z]+', '', m[1] )
        id = pepd[sequence]
        charge = spl[1]
        tuple = ( (m[0], id, charge) )
        res.append( tuple  )
    except: 
        notfound_dict[sequence] = ''
        if sequence.count('K') + sequence.count('R') < 2: print m; other +=1
        else: nontr += 1

#we have around 9511 nontryptic and 2321 other transitions
#we found 323123A => 96.5 %

q = "INSERT INTO hroest.MRMPepLink (mrm_key, peptide_key, charge)" + \
     " VALUES (%s,%s,%s) "
c.executemany( q, res)





exp_key = 3131  #yeast (all)
#exp_key = 3352  #yeast, 1200 peptides
#exp_key = 3445  #mouse (all)

#mrm atlas
all_peptide_query  = """
select distinct peptide.sequence, molecular_weight, ssrcalc,
genome_occurence, 
peptide.id as peptide_key
from peptide 
inner join peptideOrganism on peptide.id = peptideOrganism.peptide_key
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence =
peptide.sequence
where experiment_key = %s
and length( peptide.sequence ) > 1
and peptide.id in (select peptide_key from hroest.MRMPepLink)
""" % exp_key



c2 = db.cursor()
###################################
# A) store the transitions
###################################
#read in the data from DB and put into bins
# 75 peptides/sec
insert_db = True
modify_cysteins = True
read_from_db = False
c.execute( 'use ddb;' )
t = utils.db_table( c2 )
t.read( all_peptide_query )
print "executed the query"
rows = t.c.fetchall()

#1000 entries / 9 s with fast (executemany), 10.5 with index
#1000 entries / 31 s with normal (execute), 34 with index
transition_table = 'hroest.srmTransitions_yeast_mrmatlas'
peptide_table = 'hroest.srmPeptides_yeast_mrmatlas'
mass_bins = [ []  for i in range(0, 10000) ]
rt_bins = [ []  for i in range(-100, 500) ]
tmp_c  = db.cursor()
tmp2_c  = db.cursor()
progressm = progress.ProgressMeter(total=len(rows), unit='peptides')
start = time.time()
R = silver.Residues.Residues('mono')
exclude = 0
for i,row in enumerate(rows):
    progressm.update(1)
    ttmp = utils.db_table( tmp_c )
    ttmp.read( "select * from hroest.MRMAtlas_yeast_qqq_201002_q2 a inner join hroest.MRMPepLink p on p.mrm_key = a.id where peptide_key = %s order by Intensity" % row[4])
    transitions = ttmp.fetchall_groupBy( 'parent_charge')
    for ch, t_charge in transitions.iteritems(): 
        mod_sequence = ttmp.row( t_charge[0] , 'sequence').split('/')[0]
        #
        peptide = DDB.Peptide()
        peptide.set_sequence( mod_sequence )
        peptide.charge = ttmp.row( t_charge[0] , 'charge')
        peptide.ssr_calc         = t.row( row, 'ssrcalc' )
        peptide.id               = t.row( row, 'peptide_key' )
        if '[' in peptide.get_modified_sequence('SEQUEST'): continue
        #
        assert ch == peptide.charge
        S = silver.Spectrum.Spectrum(SEQUEST_mode =1 )
        S.ion_charge = peptide.charge
        S.construct_from_peptide( 
            peptide.get_modified_sequence('SEQUEST'), R.residues, 
            R.res_pairs)
        S.ass_peptide = peptide
        collider.insert_peptide_in_db(S, db, peptide_table, transition_group=i)
        prepare  = []
        vals = "type, fragment_number, group_id, q3_charge, q3 "
        q = "insert into %s (%s)" % (transition_table, vals)  + " VALUES (%s,%s,%s,%s,%s)" 
        #attention: we write the parent ion charge into the q3_charge field!
        #this to be able to still figure out the relative order of the
        #transitions later on
        for j,tt in enumerate(t_charge):
            prepare.append( [ttmp.row( tt, 'Ion_description').split('/')[0][:8],
                             j+1, i, ch, ttmp.row( tt, 'q3') ] )
        dummy = tmp2_c.executemany( q, prepare)



#important: also do the isotopes
#A2 create the additional parent ions for the isotope patterns
PATTERNS_UP_TO = 3 #up to how many isotopes should be considered
cursor = c
vals = "peptide_key, q1_charge, q1, modified_sequence, ssrcalc, isotope_nr, transition_group"
query = "SELECT parent_id, %s FROM %s" % (vals, peptide_table)
cursor.execute( query )
allpeptides =  cursor.fetchall()
prepared = []
for p in allpeptides:
    q1_charge = p[2]
    for i in range(1,1+PATTERNS_UP_TO):
        prepared.append( [p[1], p[2], p[3] + (R.mass_diffC13 * i* 1.0) / q1_charge,
                              p[4], p[5], i, p[7]] )

q = "INSERT INTO %s (%s)" % (peptide_table, vals)  + \
     " VALUES (" + "%s," *6 + "%s)"
c.executemany( q, prepared)


#END creating the MRM Atlas tables }}}
#################################################################

#################################################################
#START creating the Peptide Atlas tables {{{
from Bio import SeqIO
f = open('APD_Sc_all.fasta')
peptides = []
for p in SeqIO.parse(f, "fasta"):
    peptides.append(p.seq.tostring())

f.close()

import perl
perl.eval( "use lib '%s'" % perl_ddb_location)
import ddb #needs the perl api to get to the 2DDB
for seq in peptides:
    p = ddb.Peptide(experiment_key = 3452)
    p.peptide = seq 
    p.peptide_type = 'bioinformatics'
    p.addignore_setid()


exp_key = 3131  #yeast (all)
#pep atlas
all_peptide_query  = """
select distinct peptide.sequence, molecular_weight, ssrcalc,
genome_occurence, 
peptide.id as peptide_key
from peptide 
inner join peptideOrganism on peptide.id = peptideOrganism.peptide_key
inner join compep.ssrcalc_prediction on ssrcalc_prediction.sequence =
peptide.sequence
where experiment_key = %s
and length( peptide.sequence ) > 1
and peptide.sequence in (select sequence from ddb.peptide where experiment_key = 3452)
""" % exp_key



c2 = db.cursor()
###################################
# A) store the transitions
###################################
#read in the data from DB and put into bins
# 75 peptides/sec
insert_db = True
modify_cysteins = True
read_from_db = False
c.execute( 'use ddb;' )
t = utils.db_table( c2 )
t.read( all_peptide_query )
print "executed the query"
rows = t.c.fetchall()

#1000 entries / 9 s with fast (executemany), 10.5 with index
#1000 entries / 31 s with normal (execute), 34 with index
transition_table = 'hroest.srmTransitions_yeast_pepatlas'
peptide_table = 'hroest.srmPeptides_yeast_pepatlas'
c.execute('truncate table ' + peptide_table)
c.execute('truncate table ' + transition_table)
mass_bins = [ []  for i in range(0, 10000) ]
rt_bins = [ []  for i in range(-100, 500) ]
tmp_c  = db.cursor()
tmp2_c  = db.cursor()
progressm = progress.ProgressMeter(total=len(rows), unit='peptides')
start = time.time()
R = silver.Residues.Residues('mono')
exclude = 0
for i,row in enumerate(rows):
    progressm.update(1)
    ttmp = utils.db_table( tmp_c )
    ttmp.read( "select * from hroest.MRMAtlas_yeast_qqq_201002_q2 a inner join hroest.MRMPepLink p on p.mrm_key = a.id where peptide_key = %s order by Intensity" % row[4])
    transitions = ttmp.fetchall_groupBy( 'parent_charge')
    for ch, t_charge in transitions.iteritems(): 
        mod_sequence = ttmp.row( t_charge[0] , 'sequence').split('/')[0]
        #
        peptide = DDB.Peptide()
        peptide.set_sequence( mod_sequence )
        peptide.charge = ttmp.row( t_charge[0] , 'charge')
        peptide.ssr_calc         = t.row( row, 'ssrcalc' )
        peptide.id               = t.row( row, 'peptide_key' )
        if '[' in peptide.get_modified_sequence('SEQUEST'): continue
        #
        assert ch == peptide.charge
        S = silver.Spectrum.Spectrum(SEQUEST_mode =1 )
        S.ion_charge = peptide.charge
        S.construct_from_peptide( 
            peptide.get_modified_sequence('SEQUEST'), R.residues, 
            R.res_pairs)
        S.ass_peptide = peptide
        collider.insert_peptide_in_db(S, db, peptide_table, transition_group=i)
        prepare  = []
        vals = "type, fragment_number, group_id, q3_charge, q3 "
        q = "insert into %s (%s)" % (transition_table, vals)  + " VALUES (%s,%s,%s,%s,%s)" 
        #attention: we write the parent ion charge into the q3_charge field!
        #this to be able to still figure out the relative order of the
        #transitions later on
        for j,tt in enumerate(t_charge):
            prepare.append( [ttmp.row( tt, 'Ion_description').split('/')[0][:8],
                             j+1, i, ch, ttmp.row( tt, 'q3') ] )
        dummy = tmp2_c.executemany( q, prepare)




#END creating the Peptide Atlas tables }}}
#################################################################

#################################################################
#START creating the top3 tables {{{
observations, pI, precursor m/z, average retention time, precursor charge, peptide, protein. All peptides are proteotypic.

drop table hroest.hl_top3  ;
create table hroest.hl_top3  (
id int auto_increment primary key, 
observations int, 
total_obs int, 
pI double, 
q1 double, 
rt double, 
q1_charge int, 
sequence varchar(255),
sequence_ch varchar(255),
sequence_naked varchar(255),
protein varchar(255)
)

values = "observations , total_obs , pI , q1 , rt , q1_charge , sequence_ch, sequence , sequence_naked , protein "

import csv
r = csv.reader( open('/home/hroest/data/targets_3pep.csv'))
cursor = c
prepare  = []
for row in r:
    seq = row[5]
    mystrip = "".join( [s for s in seq if s in string.uppercase] )
    obs = row[0].split('/')[0]
    total = row[0].split('/')[1]
    rr = [obs, total]
    rr.extend( row[1:6])
    rr.append(seq.split('/')[0])
    rr.append( mystrip )
    rr.append( row[-1])
    prepare.append( rr )

c.executemany( 
" insert into hroest.hl_top3 (%s) values " % values + 
    ' (%s,%s,%s,%s,%s,%s,%s,%s,%s,%s)', prepare)


alter table hroest.hl_top3 add index(sequence_naked);
alter table hroest.hl_top3 add index(sequence);
alter table hroest.hl_top3 add index(sequence_ch);

select count(*) from hroest.hl_top3;


#946 are not tryptic
select count(*) from hroest.hl_top3 t inner join ddb.peptide p
on p.sequence = t.sequence_naked where experiment_key = 3131;

#I found 4560
select count(*) from hroest.hl_top3 t inner join hroest.srmPeptides_yeast s
on s.modified_sequence = t.sequence;

#I found 2400
select *
from hroest.hl_top3 t inner join 
hroest.MRMAtlas_yeast_qqq_201002_q2  
a on a.sequence = t.sequence_ch
group by a.sequence


#I found 4560
#create the top3 tables, which only contain all b,y ions 
#for the peptides seen in henry lams top3 estimate
create table hroest.srmPeptides_yeast_top3 as
select 
s.parent_id        ,
s.peptide_key      ,
s.modified_sequence,
s.q1_charge        ,
s.q1               ,
s.ssrcalc          ,
s.isotope_nr       ,
s.transition_group 
from hroest.hl_top3 t inner join hroest.srmPeptides_yeast s
on s.modified_sequence = t.sequence;

alter table hroest.srmPeptides_yeast_top3 add index(peptide_key);
alter table hroest.srmPeptides_yeast_top3 add index(q1);
alter table hroest.srmPeptides_yeast_top3 add index(ssrcalc);
alter table hroest.srmPeptides_yeast_top3 add index(transition_group);

create table hroest.srmTransitions_yeast_top3 as
select
srm_id          ,
group_id        ,
q3_charge       ,
q3              ,
type            ,
fragment_number 
from hroest.srmTransitions_yeast s
where s.group_id  in (select transition_group from hroest.srmPeptides_yeast_top3);
alter table hroest.srmTransitions_yeast_top3 add index(group_id);
alter table hroest.srmTransitions_yeast_top3 add index(q3);

#print "\n".join(chr(i) + " ... : Text" for i in range(65,91,2))

#END creating the top3 tables }}}
#################################################################

#################################################################
#START creating the MRM Atlas (all transitions) tables {{{
drop table hroest.srmPeptides_yeast_mrmatlasall ;
create table hroest.srmPeptides_yeast_mrmatlasall as
select  distinct
s.parent_id        ,
s.peptide_key      ,
s.modified_sequence,
s.q1_charge        ,
s.q1               ,
s.ssrcalc          ,
s.isotope_nr       ,
s.transition_group 
from  hroest.MRMPepLink t inner join hroest.srmPeptides_yeast s
on s.peptide_key = t.peptide_key;
alter table hroest.srmPeptides_yeast_mrmatlasall add index(peptide_key);
alter table hroest.srmPeptides_yeast_mrmatlasall add index(q1);
alter table hroest.srmPeptides_yeast_mrmatlasall add index(ssrcalc);
alter table hroest.srmPeptides_yeast_mrmatlasall add index(transition_group);

drop table hroest.srmTransitions_yeast_mrmatlasall ;
create table hroest.srmTransitions_yeast_mrmatlasall as
select
srm_id          ,
group_id        ,
q3_charge       ,
q3              ,
type            ,
fragment_number 
from hroest.srmTransitions_yeast s
where s.group_id  in (select transition_group from hroest.srmPeptides_yeast_mrmatlasall);
alter table hroest.srmTransitions_yeast_mrmatlasall add index(group_id);
alter table hroest.srmTransitions_yeast_mrmatlasall add index(q3);
#END creating the top3 tables }}}
#################################################################

#}}}
#################################################################
