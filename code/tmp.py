#in vim: set fdm=marker

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

#################################################################
#{{{ MRM Atlas Ion Abundance analysis

c.execute( 
"""
select id, Ion_description 
from hroest.MRMAtlas_yeast_qqq_201002_q2
"""
)
desc = c.fetchall()

d  = desc[0]

prepare = []
for d in desc:
    id = d[0]
    d = d[1]
    ion = d.split('/')[0]
    #for some reasons, 806 have a [] around them
    isotope = 0
    if ion.startswith('['): ion = ion[1:]
    if ion.find('i') != -1:
        if ion.find('i') != len(ion)-1: break
    if ion.endswith('i'): ion = ion[:-1]; isotope = 1
    charge = ion.split('^')
    if len(charge) == 1: charge = 1
    else: ion = charge[0]; charge = int( charge[1] )
    loss = ion.split('-')
    if len(loss) == 1: loss = 0
    else: ion = loss[0]; loss = int( loss[1] )
    gain = ion.split('+')
    if len(gain) == 1: gain = 0
    else: ion = gain[0]; gain = int( gain[1] )
    fragment = ion[0]
    if len(ion) > 1: nr = int(ion[1:])
    else: nr = 0
    #print  d, fragment, nr, charge, loss, gain, isotope
    prepare.append( [  id, fragment, nr, charge, loss, gain, isotope] )


create table hroest.result_mrmfragment_analysis (
id int, 
fragment char(1), 
nr int,
charge int,
loss int, 
gain int,
isotope int
)

c.executemany( 
    'insert into hroest.result_mrmfragment_analysis values (%s,%s,%s,%s,%s,%s,%s)', 
              prepare)


desc hroest.MRMAtlas_yeast_qqq_201002_q2

select sequence,Ion_description, 
@ion := LEFT(Ion_description,  INSTR( Ion_description, '/') - 1) as ion,
@fr := LEFT( @ion, 1)  as fr,
@str := LEFT(@ion,  INSTR( @ion, '-') - 1) as str
from hroest.MRMAtlas_yeast_qqq_201002_q2
limit 100;



desc hroest.result_mrmfragment_analysis;
#there is no fragment a with any gains/losses
select @all := count(*)  from hroest.result_mrmfragment_analysis;
select count(*) / @all * 100 , fragment from hroest.result_mrmfragment_analysis
group by fragment
;
select count(*) / @all * 100, loss from hroest.result_mrmfragment_analysis
group by loss
;
select count(*) / @all * 100, gain from hroest.result_mrmfragment_analysis
group by gain
;
drop table if exists tmp ;
create temporary table tmp as
select count(*) / @all * 100 as pcnt, fragment, nr from hroest.result_mrmfragment_analysis
group by fragment, nr
;
select  * from tmp where pcnt > 1.0;
#
drop table if exists tmp ;
create temporary table tmp as
select count(*) / @all * 100 as pcnt, fragment, nr, loss, gain 
from hroest.result_mrmfragment_analysis
group by fragment, nr, loss, gain
;
select  * from tmp where pcnt > 0.1;

#interestingly, a loss of 17 and 18 is nearly twice as likely from y than b
select count(*) / @all * 100, fragment, loss from hroest.result_mrmfragment_analysis
group by fragment, loss
;
select count(*) / @all * 100, fragment, charge from hroest.result_mrmfragment_analysis
group by fragment, charge
;

select count(*) / @all * 100, fragment, gain from hroest.result_mrmfragment_analysis
group by fragment, gain
;

#how often do the different charge combinations happen?
select count(*) / @all * 100, charge, parent_charge 
from hroest.result_mrmfragment_analysis an
inner join  MRMAtlas_yeast_qqq_201002_q2 at on at.id = an.id
group by charge, parent_charge
;


select @Intensity := 1000;
select @all := count(*) from MRMAtlas_yeast_qqq_201002_q2
where Intensity > @Intensity; 
select count(*) / @all * 100 from hroest.result_mrmfragment_analysis an
inner join  MRMAtlas_yeast_qqq_201002_q2 at on at.id = an.id
where fragment in ('b', 'y')
and loss = 0 and gain = 0
and charge in (1)
and isotope  = 0
and Intensity > @Intensity; 
;

select count(*) / @all * 100 from hroest.result_mrmfragment_analysis
where 
isotope  = 1
;

+-----------------------+----------+
| count(*) / @all * 100 | fragment |
+-----------------------+----------+
|                0.2141 | ?        | 
|                4.1522 | a        | 
|               15.1250 | b        | 
|                0.5977 | p        | 
|               79.9110 | y        | 
+-----------------------+----------+

-----------------------+------+
| count(*) / @all * 100 | loss |
+-----------------------+------+
|               90.3802 |    0 | 
|                3.3425 |   17 | 
|                3.5521 |   18 | 
|                0.0427 |   34 | 
|                0.8073 |   35 | 
|                0.5550 |   36 | 
|                0.5881 |   44 | 
|                0.1346 |   45 | 
|                0.5768 |   46 | 
|                0.0069 |   64 | 
|                0.0003 |   82 | 
|                0.0134 |   91 | 
+-----------------------+------+

+-----------------------+------+
| count(*) / @all * 100 | gain |
+-----------------------+------+
|               99.6241 |    0 | 
|                0.3759 |   18 | 
+-----------------------+------+


+---------+----------+------+
| pcnt    | fragment | nr   |
+---------+----------+------+
|  2.3355 | b        |    4 | 
|  3.2401 | b        |    5 | 
|  2.6616 | b        |    6 | 
|  2.0773 | b        |    7 | 
|  1.6513 | b        |    8 | 
|  1.1225 | b        |    9 | 
|  1.9779 | y        |    3 | 
|  8.9266 | y        |    4 | 
| 10.6080 | y        |    5 | 
| 10.9603 | y        |    6 | 
| 11.2863 | y        |    7 | 
| 10.6573 | y        |    8 | 
|  8.7869 | y        |    9 | 
|  6.6570 | y        |   10 | 
|  4.6102 | y        |   11 | 
|  2.6183 | y        |   12 | 
|  1.2190 | y        |   13 | 
+---------+----------+------+


+-----------------------+----------+--------+
| count(*) / @all * 100 | fragment | charge |
+-----------------------+----------+--------+
|                0.2141 | ?        |      1 | 
|                2.9073 | a        |      1 | 
|                1.1055 | a        |      2 | 
|                0.1394 | a        |      3 | 
|               13.1459 | b        |      1 | 
|                1.7713 | b        |      2 | 
|                0.2078 | b        |      3 | 
|                0.0009 | p        |      1 | 
|                0.5789 | p        |      2 | 
|                0.0179 | p        |      3 | 
|               69.3401 | y        |      1 | 
|                9.9658 | y        |      2 | 
|                0.6052 | y        |      3 | 
+-----------------------+----------+--------+


select @all / count(distinct sequence) from MRMAtlas_yeast_qqq_201002_q2;
select @all / count(distinct peptide_key) from MRMAtlas_yeast_qqq_201002_q2 a
inner join MRMPepLink l on l.mrm_key = a.id;

#the average minimal intensity is 1704
select sum(minint) / 44215 from 
( select min(Intensity) as minint, sequence from MRMAtlas_yeast_qqq_201002_q2
 group by sequence
) tmp
;

#}}}
{{{#some analysis of the peptide atlas data

#peptide Atlas
#has key 3452
cursor.execute(
"""
create table hroest.delete as select *  from ddb.peptide where experiment_key = 3452;
create table hroest.patlasNoGenome as select * from hroest.delete d where
    d.sequence not in (select sequence from ddb.peptide where experiment_key =
    3131);
"""
)
cursor.execute("select sequence from hroest.patlasNoGenome")
list = cursor.fetchall()
cursor.execute("select sequence from hroest.patlasNoGenome")

nontr = 0
other = []
for s in list:
    seq = s[0]
    if seq.count('K') + seq.count('R') > 1: nontr += 1
    else: other.append( seq)
    

#total 25075
select count(*) from hroest.patlasNoGenome; 
#22139 are tryptic
select count(*) from hroest.patlasNoGenome where
sequence like '%K' or
sequence like '%R' ;
#6540 contain no missed cleavages at the end
drop table hroest.patlasNoGenome2 ;
create table hroest.patlasNoGenome2 as select * 
#select count(*) 
from hroest.patlasNoGenome where
(sequence like '%K' or
sequence like '%R' ) and not  (
sequence like '%K%K' or
sequence like '%K%R' or
sequence like '%R%R' or
sequence like '%R%K');

cursor.execute("select sequence from hroest.patlasNoGenome2;")
nf = cursor.fetchall()
cursor.execute("select sequence from ddb.peptide where experiment_key = 3131")
all = cursor.fetchall()

count = 0
others = []
for n in nf:
    ishere = False
    for a in all:
        if a[0].find(n[0]) != -1: print a, n; ishere = True
    if ishere: count += 1
    else: others.append( n[0] )


len(nf)
count

DIENQYETQITQIEHEVSSSGQEVQSSAK
cursor.execute("drop table if exists hroest.patlasNoGenome3" )
cursor.execute("create table hroest.patlasNoGenome3 (sequence varchar(255) )" )
cursor.executemany( 'insert into hroest.patlasNoGenome3 (sequence) values (%s)', others)

#225 are human contaminants
#select count(*) from ddb.peptide 
drop table hroest.patlasNoGenome4 ;
create table hroest.patlasNoGenome4 as
select sequence from hroest.patlasNoGenome3 pa
where pa.sequence not in (select sequence from ddb.peptide where experiment_key = 3130)




select * from hroest.patlasNoGenome limit 300;






}}}
#################################################################

#################################################################
{{{ #Analysis of the Manndata

#make map table
drop table if exists tmp_pepmapping;
create temporary table tmp_pepmapping as
select distinct peptide_key, modified_sequence from 
hroest.srmPeptides_yeast ;
alter table tmp_pepmapping add index(peptide_key);
alter table tmp_pepmapping add index(modified_sequence);


#desc  manndata_lfq;
#select count(*) from manndata_lfq;
drop table if exists tmp_peptide ; 
create temporary table tmp_peptide as 
select peptide, sum(obs) as nr, max(average) as av 
from manndata_lfq group by peptide;





#select * from tmp_peptide limit 100;

#in total 713'715 features mapped to identifications
select sum(obs) from manndata_lfq;
#we have 33k peptides
select count(*) from tmp_peptide;
#we have 26k tryptic peptides
select count(*) from tmp_peptide t 
inner join hroest.tmp_pepmapping p on
t.peptide = p.modified_sequence ;
#24k proteotypic peptides
select count(*) from tmp_peptide t 
inner join hroest.tmp_pepmapping p on t.peptide = p.modified_sequence 
inner join ddb.peptideOrganism o on o.peptide_key = p.peptide_key
where genome_occurence = 1
;


#we have 3687 unique proteins
#with maxquant they have 

drop table tmp_allfoundproteins ;
create temporary table tmp_allfoundproteins as 
select count(*) as nr_peptides, protein_key from  tmp_peptide t 
inner join tmp_pepmapping map on map.modified_sequence = t.peptide 
inner join ddb.peptide p on map.peptide_key = p.id
inner join ddb.peptideOrganism o on o.peptide_key = p.id
inner join ddb.protPepLink l on l.peptide_key = p.id
where experiment_key = 3131
and genome_occurence = 1
group by protein_key
order by count(*);

#we have 2305 proteins of which we have 3 or more tryptic peptides
select count(*) from tmp_allfoundproteins where nr_peptides > 2;
select count(*) from tmp_allfoundproteins where nr_peptides = 2;
select count(*) from tmp_allfoundproteins where nr_peptides = 1;

#estimate FDR
select count(*) *0.01*0.01*0.01 from tmp_allfoundproteins where nr_peptides > 2;
select count(*) *0.01*0.01 from tmp_allfoundproteins where nr_peptides = 2;
select count(*) *0.01 from tmp_allfoundproteins where nr_peptides = 1;

#from 27k peptides, we have 8k in the lab
drop table if exists tmp_mannanalysis_pepperprot ;
create table tmp_mannanalysis_pepperprot as
select t.peptide, av as intensity, genome_occurence, protein_key from  tmp_peptide t 
inner join tmp_pepmapping map on map.modified_sequence = t.peptide 
inner join ddb.peptide p on map.peptide_key = p.id
inner join ddb.peptideOrganism o on o.peptide_key = p.id
inner join ddb.protPepLink l on l.peptide_key = p.id
#inner join lab_peptides lab on p.sequence = lab.sequence
where experiment_key = 3131
and genome_occurence = 1
order by protein_key, intensity
;




cursor.execute( 'select * from hroest.tmp_mannanalysis_pepperprot')
mypep = cursor.fetchall()

cursor.execute( 'select * from hroest.lab_peptides')
lab = cursor.fetchall()

labs = [l[3] for l in lab]

proteins = {}
for p in mypep:
    if proteins.has_key(p[3]): proteins[p[3]].append( p)
    else: proteins[p[3]] = [p]


for nr in range(1,7):
    intop = []
    for k,v in proteins.iteritems():
        intop.extend( v[-nr:])
    #we have 9234 of the top3 peptides
    #679 + 2*469 + 2539*3
    haveinlab = 0
    for p in intop:
        if p[0] in labs: haveinlab += 1
    print nr, haveinlab, len(intop), haveinlab*1.0 / len(intop)

#of those 9234 we have 4393

#of the top1 (3687) we have 1903 = 51%
#of the top2 (6695) we have 3326 = 50%
#of the top3 (9234) we have 4393 = 47%
#of the top4 (11420) we have 5186 = 45%
#of the top5 (13313) we have 5821 = 43%
#of the top6 (14965) we have 5276 = 41%


#TODO run mayo on my data:
    go to peptides.csv, get all proteotypic peptides (only one n_proteins)
    then write a file with 
    specid, peptide, protid, modification, score
    dummy,  PEPTIDE, DECOY_da, none, dummy

=> then estimate FDR of the protein identifications





#tinas proteins
#select sequence_key, ac from ddb.isbAc where ac in ( 
select ac, t.peptide, 
t.av as max_triplicate_log_intensity,
t.nr as nr_files_observed, 
genome_occurence
from  tmp_peptide t 
inner join tmp_pepmapping map on map.modified_sequence = t.peptide 
inner join ddb.peptide p on map.peptide_key = p.id
inner join ddb.peptideOrganism o on o.peptide_key = p.id
inner join ddb.protPepLink l on l.peptide_key = p.id
inner join ddb.protein prot on l.protein_key = prot.id
inner join ddb.isbAc ac on ac.sequence_key = prot.sequence_key
where p.experiment_key = 3131
and ac.ac in 
(
'YBR249C',
'YCL017C',
'YDR283C',
'YEL031W',
'YER178W',
'YFL030W',
'YGL006W',
'YGL248W',
'YGL253W',
'YGR193C',
'YGR240C',
'YHR107C',
'YHR166C',
'YHR183W',
'YIL084C',
'YJL026W',
'YJL136C',
'YJR051W',
'YKL060C',
'YKL145W',
'YKR031C',
'YLR039C',
'YLR058C',
'YLR249W',
'YLR403W',
'YML109W',
'YMR037C',
'YMR170C',
'YMR205C',
'YMR220W',
'YNR067C',
'YOL116W',
'YPR118W')
group by ac,peptide
order by ac, av DESC
;



#now check the peptide atlas
select p.sequence_key, ac, pep.sequence from ddb.isbAc  ac
inner join ddb.protein p  on p.sequence_key = ac.sequence_key
inner join ddb.protPepLink l on l.protein_key = p.id
inner join ddb.peptide pep on l.peptide_key = pep.id
inner join ddb.peptideOrganism o on o.peptide_key = pep.id
where p.experiment_key = 3131
and genome_occurence = 1
and  pep.sequence in (select sequence from ddb.peptide where experiment_key = 3452)
and ac in ( 
'YBR249C',
'YCL017C',
'YDR283C',
'YEL031W',
'YER178W',
'YFL030W',
'YGL006W',
'YGL248W',
'YGL253W',
'YGR193C',
'YGR240C',
'YHR107C',
'YHR166C',
'YHR183W',
'YIL084C',
'YJL026W',
'YJL136C',
'YJR051W',
'YKL060C',
'YKL145W',
'YKR031C',
'YLR039C',
'YLR058C',
'YLR249W',
'YLR403W',
'YML109W',
'YMR037C',
'YMR170C',
'YMR205C',
'YMR220W',
'YNR067C',
'YOL116W',
'YPR118W')
group by ac, pep.sequence
order by ac
;




KPADLASLLLNSAGDAQGDEAPALK                   
| YEL031W | LANVSAVTNIIR                    
| YEL031W | FQFSSALK                        
| YEL031W | GTVITPEIR                       
| YEL031W | SDDNQLLFR                       



select average, stdev, obs, file from manndata_lfq 
where peptide = 'SDDNQLLFR'
and average != 0
and obs = 3
;


select average, stdev, obs, file from manndata_lfq 
where peptide = 'GTVITPEIR'
and average != 0
#and obs = 3
;


}}}
#################################################################



###########################################################################
### some code to calculate where we were without RT information {{{
##NO rt information
#

query = """
select q1, q3 
from %(peptable)s
inner join %(tratable)s on parent_id = parent_key
where 
q1 between 400 and 1200 
and q3 between 400 and 1200 
""" % {'peptable' : par.peptide_table ,
       'tratable' : par.transition_table }
cursor.execute( query)

lines = cursor.fetchall()
mybins = [ [0 for i in range(0, 1200) ] for i in range(0, 1200 / 0.7+ 1 ) ]
for l in lines:
    q1 = int( numpy.floor( l[0]/ 0.7 ))
    q3 = int( numpy.floor( l[1] ))
    mybins[q1][q3] += 1


total = 0
nr_bins = 0
mymax = 0
single_spots = 0
mybest = (0,0)
for out in range(400, 1200/0.7+1):
    for inn in range(400, 1200):
        total += mybins[out][inn]
        if mybins[out][inn] > mymax:
            mymax = mybins[out][inn]
            mybest = (out, inn)
        nr_bins += 1
        if mybins[out][inn] == 1: single_spots += 1

avg = 1.0* total / nr_bins
unique_pct = single_spots * 1.0 / total
}}}


